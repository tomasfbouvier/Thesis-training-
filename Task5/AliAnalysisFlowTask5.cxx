/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

 /**************************************************************************
  * Example analysis task for flow calculations                            *
  *                                                                        *
  * Author: Zuzana Moravcova, NBI, 2019 --                                 *
  *         zuzana.moravcova@cern.ch                                       *
  **************************************************************************/

/* AliAnalysisFlowTask5
 *
 * example flow class:
 * event selection
 * track selection
 * loading particle weights if needed (possible to change run (in)dependent)
 * filling Q vectors (only needed combinations of harmonics & powers! check if this is what you need!)
 * calculate m-particle correlations (using generic algorithm)
 */

#include "AliAnalysisFlowTask5.h"

class AliAnalysisFlowTask5;

using namespace std;

ClassImp(AliAnalysisFlowTask5)

AliAnalysisFlowTask5::AliAnalysisFlowTask5() : AliAnalysisTaskSE(),
    fAOD(0),
    fTrigger(AliVEvent::kINT7),
    fIsPP(kFALSE),
    fEventRejectAddPileUp(kFALSE),
    fFilterBit(96),
    fCentBinNum(10),
    fCentMin(0),
    fCentMax(100),
    fCentrality(0.0),
    fPtMin(0.2),
    fPtMax(5.0),
    fAbsEtaMax(0.8),
    fPVtxCutZ(10.0),
    fPVz(99.9),
    fCentEstimator("V0M"),
    fOutputList(0),
    fhEventCounter(0),
    fDCAtoVtxZ(3.0),
    fDCAtoVtxXY(3.0),
    fNTPCClusters(70),
    fNITSClusters(2),
    fTPCchi2perClusterMin(0.1),
    fTPCchi2perClusterMax(4.)
    {}
//_____________________________________________________________________________
AliAnalysisFlowTask5::AliAnalysisFlowTask5(const char* name, Bool_t bUseWeights) : AliAnalysisTaskSE(name),
    fAOD(0),
    fTrigger(AliVEvent::kINT7),
    fIsPP(kFALSE),
    fEventRejectAddPileUp(kFALSE),
    fUseWeights(bUseWeights),
    fFilterBit(96),
    fCentBinNum(10),
    fCentMin(0),
    fCentMax(100),
    fCentrality(0.0),
    fPtMin(0.2),
    fPtMax(5.0),
    fAbsEtaMax(0.8),
    fPVtxCutZ(10.0),
    fPVz(99.9),
    fCentEstimator("V0M"),
    fOutputList(0),
    fInputWeights(0),
    fhEventCounter(0),
    fhWeight(0),
    fgap(1.),
    fDCAtoVtxZ(3.0),
    fDCAtoVtxXY(3.0),
    fNTPCClusters(70),
    fNITSClusters(2),
    fTPCchi2perClusterMin(0.1),
    fTPCchi2perClusterMax(4.),
    fGapping(kFALSE),
    fWeighting(kTRUE),
    fDiff(kFALSE),
    maxHarm(5),
    maxWeightPower(5),
    fFlowVecQpos{},
    fFlowVecQneg{}
{
    DefineInput(0, TChain::Class());
    if(bUseWeights) { DefineInput(1, TList::Class()); }
    DefineOutput(1, TList::Class());
}
//_____________________________________________________________________________
AliAnalysisFlowTask5::~AliAnalysisFlowTask5()
{
    if(fOutputList) delete fOutputList;
    if(fhWeight) delete fhWeight;
}
//_____________________________________________________________________________
void AliAnalysisFlowTask5::UserCreateOutputObjects()
{

    // OpenFile(1);

    // creating output lists
    fOutputList = new TList();
    fOutputList->SetOwner(kTRUE);
    
    
    
    TString sEventCounterLabel[] = {"Input", "Trigger OK", "Event cuts OK", "Event OK"};
    const Int_t iEventCounterBins = sizeof(sEventCounterLabel)/sizeof(sEventCounterLabel[0]);
    fhEventCounter = new TH1D("fhEventCounter","Event Counter",iEventCounterBins,0,iEventCounterBins);
    fOutputList->Add(fhEventCounter);

    const Double_t bin_edges[] = {0.,5.,10.,20.,30.,40.,50.};
    
    std::string names2[]= {"_1","_2","_3","_4","_5","_6","_7","_8","_9", "_10", "MEAN"};
    
    for(Int_t i(0); i<11; i++){
        fProfiles_collection_SC_5_2[i]= new TProfile(("4pc_5_2"+ names2[i]).c_str(),"Profile eta vs value; centrality;  <4>(5,2) ",6, bin_edges);
        fProfiles_collection_SC_5_3[i]= new TProfile(("4pc_5_3"+ names2[i]).c_str(),"Profile eta vs value; centrality;  <4>(5,3) ",6, bin_edges);
        fProfiles_collection_SC_4_3[i]= new TProfile(("4pc_4_3"+ names2[i]).c_str(),"Profile eta vs value; centrality;  <4>(4,3) ",6, bin_edges);
        fProfiles_collection_SC_4_2[i]= new TProfile(("4pc_4_2"+ names2[i]).c_str(),"Profile eta vs value; centrality;  <4>(4,2) ",6, bin_edges);
        fProfiles_collection_SC_3_2[i]= new TProfile(("4pc_3_2"+ names2[i]).c_str(),"Profile eta vs value; centrality;  <4>(3,2) ",6, bin_edges);
        
        fProfiles_collection_2pc_2[i] = new TProfile(("2pc_2"+ names2[i]).c_str(),"Profile eta vs value; centrality;  <2>(2) ",6, bin_edges);
        fProfiles_collection_2pc_3[i] = new TProfile(("2pc_3"+ names2[i]).c_str(),"Profile eta vs value; centrality;  <2>(2) ",6, bin_edges);
        fProfiles_collection_2pc_4[i] = new TProfile(("2pc_4"+ names2[i]).c_str(),"Profile eta vs value; centrality;  <2>(2) ",6, bin_edges);
        fProfiles_collection_2pc_5[i] = new TProfile(("2pc_5"+ names2[i]).c_str(),"Profile eta vs value; centrality;  <2>(2) ",6, bin_edges);
        
        fOutputList->Add(fProfiles_collection_SC_5_2[i]);
        fOutputList->Add(fProfiles_collection_SC_5_3[i]);
        fOutputList->Add(fProfiles_collection_SC_4_3[i]);
        fOutputList->Add(fProfiles_collection_SC_4_2[i]);
        fOutputList->Add(fProfiles_collection_SC_3_2[i]);
        
        fOutputList->Add(fProfiles_collection_2pc_2[i]);
        fOutputList->Add(fProfiles_collection_2pc_3[i]);
        fOutputList->Add(fProfiles_collection_2pc_4[i]);
        fOutputList->Add(fProfiles_collection_2pc_5[i]);
        
    }

    if(fUseWeights) fInputWeights = (TList*) GetInputData(1);
    
    PostData(1, fOutputList);
}
//_____________________________________________________________________________
void AliAnalysisFlowTask5::UserExec(Option_t *)
{
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if(!fAOD) { return; }

    fhEventCounter->Fill("Input",1);

    if(!IsEventSelected()) { return; }

    fhEventCounter->Fill("Event OK",1);
    
    if(fUseWeights && !AreWeightsLoaded()) {
        return; }
    

    
    Int_t iTracks(fAOD->GetNumberOfTracks());
    if(iTracks < 1 ) { return; }
    
    ResetFlowVector(fFlowVecQpos, maxHarm, maxWeightPower);
    //!ResetFlowVector(fFlowVecQneg, maxHarm, maxWeightPower);
    
    for(Int_t i(0); i < iTracks; i++) {
        
        AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));
        if(!track || !IsTrackSelected(track)) { continue; }
        
        if(fWeighting) weight =  fhWeight->GetBinContent(fhWeight->GetBin(track->Phi()));
        else weight=1.;
        
        if(!fGapping) // no eta gap
              {
                for(Int_t iHarm(0); iHarm <= maxHarm; iHarm++){
                  for(Int_t iPower(0); iPower <= maxWeightPower; iPower++)
                  {
                    Double_t dCos = TMath::Power(weight,iPower) * TMath::Cos(iHarm * track->Phi());
                    Double_t dSin = TMath::Power(weight,iPower) * TMath::Sin(iHarm * track->Phi());
                    fFlowVecQpos[iHarm][iPower] += TComplex(dCos,dSin,kFALSE);
                  }
                }
              }
        
        else
              {
                // RFP in positive eta acceptance
                if(track->Eta() > fgap/2.)
                {
                  for(Int_t iHarm(0); iHarm <= maxHarm; iHarm++){
                    for(Int_t iPower(0); iPower <= maxWeightPower; iPower++)
                    {
                      Double_t dCos = TMath::Power(weight,iPower) * TMath::Cos(iHarm * track->Phi());
                      Double_t dSin = TMath::Power(weight,iPower) * TMath::Sin(iHarm * track->Phi());
                      fFlowVecQpos[iHarm][iPower] += TComplex(dCos,dSin,kFALSE);
                    }
                  }
                }
                // RFP in negative eta acceptance
                if(track->Eta() < - fgap/2.)
                {
                  for(Int_t iHarm(0); iHarm <= maxHarm; iHarm++){
                    for(Int_t iPower(0); iPower <= maxWeightPower; iPower++)
                    {
                      Double_t dCos = TMath::Power(weight,iPower) * TMath::Cos(iHarm * track->Phi());
                      Double_t dSin = TMath::Power(weight,iPower) * TMath::Sin(iHarm * track->Phi());
                      fFlowVecQneg[iHarm][iPower] += TComplex(dCos,dSin,kFALSE);
                    }
                  }
                }
             }
    }

    fiHarm[0]=2;fiHarm[1]=5;fiHarm[2]=-5;fiHarm[3]=-2;
    
    Calculate_correlation_4(fProfiles_collection_SC_5_2);
    
    fiHarm[0]=3;fiHarm[1]=5;fiHarm[2]=-5;fiHarm[3]=-3;
    
    Calculate_correlation_4(fProfiles_collection_SC_5_3);
    
    fiHarm[0]=3;fiHarm[1]=4;fiHarm[2]=-4;fiHarm[3]=-3;
    
    Calculate_correlation_4(fProfiles_collection_SC_4_3);
    
    fiHarm[0]=2;fiHarm[1]=4;fiHarm[2]=-4;fiHarm[3]=-2;
    
    Calculate_correlation_4(fProfiles_collection_SC_4_2);
    
    fiHarm[0]=2;fiHarm[1]=3;fiHarm[2]=-2;fiHarm[3]=-3;
    
    Calculate_correlation_4(fProfiles_collection_SC_3_2);
    
    fiHarm[0]=2;fiHarm[1]=-2;
    
    Calculate_correlation_2(fProfiles_collection_2pc_2);
    
    fiHarm[0]=3;fiHarm[1]=-3;
    
    Calculate_correlation_2(fProfiles_collection_2pc_3);
    
    fiHarm[0]=4;fiHarm[1]=-4;
    
    Calculate_correlation_2(fProfiles_collection_2pc_4);
    
    fiHarm[0]=5;fiHarm[1]=-5;
    
    Calculate_correlation_2(fProfiles_collection_2pc_5);
    
    PostData(1, fOutputList);
    
    
}
//_____________________________________________________________________________
void AliAnalysisFlowTask5::Terminate(Option_t *)
{}
//_____________________________________________________________________________
Bool_t AliAnalysisFlowTask5::IsEventSelected()
{
    AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler* inputHandler = (AliInputEventHandler*) mgr->GetInputEventHandler();
    UInt_t fSelectMask = inputHandler->IsEventSelected();
    if(!(fSelectMask & fTrigger)) { return kFALSE; }
    fhEventCounter->Fill("Trigger OK",1);

    AliMultSelection* multSelection = (AliMultSelection*) fAOD->FindListObject("MultSelection");
    if(!multSelection) { AliError("AliMultSelection object not found! Returning -1"); return -1; }
    fhEventCounter->Fill("Multiplicity OK",1);

    if(fIsPP) fEventCuts.OverrideAutomaticTriggerSelection(AliVEvent::kHighMultV0, true);
    if(!fEventCuts.AcceptEvent(fAOD)) { return kFALSE; }
    fhEventCounter->Fill("Event cuts OK",1);

    Float_t dPercentile = multSelection->GetMultiplicityPercentile(fCentEstimator);
    if(dPercentile > 100 || dPercentile < 0) { AliWarning("Centrality percentile estimated not within 0-100 range. Returning -1"); return -1; }
    if(fEventRejectAddPileUp && !fIsPP && dPercentile > 0 && dPercentile < 10 && IsEventRejectedAddPileUp()) { return kFALSE; }
    fhEventCounter->Fill("Centrality OK",1);
    fCentrality = (Double_t) dPercentile;

    fPVz = fAOD->GetPrimaryVertex()->GetZ();
    if(TMath::Abs(fPVz) >= fPVtxCutZ) { return kFALSE; }
    fhEventCounter->Fill("PVzOK",1);
    


    return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliAnalysisFlowTask5::IsTrackSelected(const AliAODTrack* track) const
{
    if(!track->TestFilterBit(fFilterBit)) { return kFALSE; }
    if(track->GetTPCNcls() < 70 && fFilterBit != 2) { return kFALSE; }
    if(fPtMin > 0 && track->Pt() < fPtMin) { return kFALSE; }
    if(fPtMax > 0 && track->Pt() > fPtMax) { return kFALSE; }
    if(fAbsEtaMax > 0 && TMath::Abs(track->Eta()) > fAbsEtaMax) { return kFALSE; }

    return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliAnalysisFlowTask5::IsEventRejectedAddPileUp() const
{
  // Check for additional pile-up rejection in Run 2 Pb-Pb collisions (15o, 17n)
  // based on multiplicity correlations
  // ***************************************************************************

  Bool_t bIs17n = kFALSE;
  Bool_t bIs15o = kFALSE;

  Int_t iRunNumber = fAOD->GetRunNumber();
  if(iRunNumber >= 244824 && iRunNumber <= 246994) { bIs15o = kTRUE; }
  else if(iRunNumber == 280235 || iRunNumber == 20234) { bIs17n = kTRUE; }
  else { return kFALSE; }

  // recounting multiplcities
  const Int_t multESD = ((AliAODHeader*) fAOD->GetHeader())->GetNumberOfESDTracks();
  const Int_t nTracks = fAOD->GetNumberOfTracks();
  Int_t multTPC32 = 0;
  Int_t multTPC128 = 0;
  Int_t multTOF = 0;
  Int_t multTrk = 0;
  Double_t multESDTPCdif = 0.0;
  Double_t v0Centr = 0.0;

  for(Int_t it(0); it < nTracks; it++)
  {
    AliAODTrack* track = (AliAODTrack*) fAOD->GetTrack(it);
    if(!track) { continue; }

    if(track->TestFilterBit(32))
    {
      multTPC32++;
      if(TMath::Abs(track->GetTOFsignalDz()) <= 10.0 && track->GetTOFsignal() >= 12000.0 && track->GetTOFsignal() <= 25000.0) { multTOF++; }
      if((TMath::Abs(track->Eta())) < fAbsEtaMax && (track->GetTPCNcls() >= 70) && (track->Pt() >= fPtMin) && (track->Pt() < fPtMax)) { multTrk++; }
    }

    if(track->TestFilterBit(128)) { multTPC128++; }
  }

  if(bIs17n)
  {
    multESDTPCdif = multESD - (6.6164 + 3.64583*multTPC128 + 0.000126397*multTPC128*multTPC128);
    if(multESDTPCdif > 1000) { return kTRUE; }
    if( ((AliAODHeader*) fAOD->GetHeader())->GetRefMultiplicityComb08() < 0) { return kTRUE; }
  }

  if(bIs15o)
  {
    multESDTPCdif = multESD - 3.38*multTPC128;
    if(multESDTPCdif > 500) { return kTRUE; }

    TF1 fMultTOFLowCut = TF1("fMultTOFLowCut", "[0]+[1]*x+[2]*x*x+[3]*x*x*x - 4.*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x+[9]*x*x*x*x*x)", 0, 10000);
    fMultTOFLowCut.SetParameters(-1.0178, 0.333132, 9.10282e-05, -1.61861e-08, 1.47848, 0.0385923, -5.06153e-05, 4.37641e-08, -1.69082e-11, 2.35085e-15);
    if(Double_t(multTOF) < fMultTOFLowCut.Eval(Double_t (multTPC32))) { return kTRUE; }

    TF1 fMultTOFHighCut = TF1("fMultTOFHighCut", "[0]+[1]*x+[2]*x*x+[3]*x*x*x + 4.*([4]+[5]*x+[6]*x*x+[7]*x*x*x+[8]*x*x*x*x+[9]*x*x*x*x*x)", 0, 10000);
    fMultTOFHighCut.SetParameters(-1.0178, 0.333132, 9.10282e-05, -1.61861e-08, 1.47848, 0.0385923, -5.06153e-05, 4.37641e-08, -1.69082e-11, 2.35085e-15);
    if(Double_t(multTOF) > fMultTOFHighCut.Eval(Double_t (multTPC32))) { return kTRUE; }

    AliMultSelection* multSelection = (AliMultSelection*) fAOD->FindListObject("MultSelection");
    if(!multSelection) { AliError("AliMultSelection object not found! Returning -1"); return -1; }
    v0Centr = multSelection->GetMultiplicityPercentile("V0M");

    TF1 fMultCentLowCut = TF1("fMultCentLowCut", "[0]+[1]*x+[2]*exp([3]-[4]*x) - 5.*([5]+[6]*exp([7]-[8]*x))", 0, 100);
    fMultCentLowCut.SetParameters(-6.15980e+02, 4.89828e+00, 4.84776e+03, -5.22988e-01, 3.04363e-02, -1.21144e+01, 2.95321e+02, -9.20062e-01, 2.17372e-02);
    if(Double_t(multTrk) < fMultCentLowCut.Eval(v0Centr)) { return kTRUE; }
  }

  return kFALSE;
}

Bool_t AliAnalysisFlowTask5::AreWeightsLoaded()
{
  if(!fInputWeights) {AliError("Weights input list not loaded"); return kFALSE; }
  fhWeight = (TH1D*)fInputWeights->FindObject("hWeight");
  if(!fhWeight) {AliError("Weights not loaded"); return kFALSE; }
  return kTRUE;
}
//_____________________________________________________________________________
Double_t AliAnalysisFlowTask5::GetWeight(const Double_t dPhi)
{
  if(!fUseWeights) return 1.0;
  if(!fhWeight) {AliFatal("Weights not loaded"); return -1.0; }
  return fhWeight->GetBinContent(fhWeight->FindFixBin(dPhi));
}

// ============================================================================

TComplex AliAnalysisFlowTask5::Q(const Int_t n, const Int_t p) const
{
 if (n < 0) return TComplex::Conjugate(fFlowVecQpos[-n][p]);
 else return fFlowVecQpos[n][p];
}
// ============================================================================
TComplex AliAnalysisFlowTask5::QGapPos(const Int_t n, const Int_t p) const
{
 if (n < 0) return TComplex::Conjugate(fFlowVecQpos[-n][p]);
 else return fFlowVecQpos[n][p];
}
// ============================================================================

void AliAnalysisFlowTask5::ResetFlowVector(TComplex (&array)[fFlowNumHarmonicsMax][fFlowNumWeightPowersMax], Int_t maxHarm, Int_t maxWeightPower)
{
  // Reset RFPs (Q) array values to TComplex(0,0,kFALSE) for given array
  // *************************************************************
  for(Int_t iHarm(0); iHarm <= maxHarm; ++iHarm) {
    for(Int_t iPower(0); iPower <= maxWeightPower; ++iPower) {
      array[iHarm][iPower](0.0,0.0);
    }
  }
  return;
}


void AliAnalysisFlowTask5::Calculate_correlation_2(TProfile*  filling_collect[]){
    

   //! cout<<cNom<<"\t";
    // important to distinguish between different "setup"

    Double_t value = 0.0;
    Double_t valueNeg = 0.0;
    
    TComplex cNom = TComplex(0.0,0.0,kFALSE);
    TComplex cDenom = TComplex(0.0,0.0,kFALSE);
    TComplex cNomNeg = TComplex(0.0,0.0,kFALSE);
    TComplex cDenomNeg = TComplex(0.0,0.0,kFALSE);
    
    if(!fGapping) { // no gap
        cDenom = Two(0,0);
        cNom = Two(fiHarm[0],fiHarm[1]);
     
      }
    else { // has gap
      cDenom = TwoGap(0,0);
      cNom = TwoGap(fiHarm[0],fiHarm[1]);
          
    }
    
    Bool_t bFillPos = kFALSE;
    Bool_t bFillNeg = kFALSE;
    
    if(cDenom.Re() > 0.0) { bFillPos = kTRUE; value = cNom.Re() / cDenom.Re(); }
    if(bFillPos && TMath::Abs(value) > 1.0) { bFillPos = kFALSE; }
    
    if(value > -100 and fCentrality< 90){
        Int_t random_number = rand() % 10 ;
        if(bFillPos) filling_collect[10]->Fill(fCentrality, value, cDenom.Re());
        if(bFillPos) filling_collect[random_number]->Fill(fCentrality, value, cDenom.Re());
        //!fMyProfile_v_SC_5_2_2->Fill(fCentrality,value, cDenom.Re());
    }
    
    return ;
    
}


void AliAnalysisFlowTask5::Calculate_correlation_4(TProfile*  filling_collect[]){
    

   //! cout<<cNom<<"\t";
    // important to distinguish between different "setup"
    
    Double_t value=0.;

    TComplex cNom = TComplex(0.0,0.0,kFALSE);
    TComplex cDenom = TComplex(0.0,0.0,kFALSE);
    TComplex cNomNeg = TComplex(0.0,0.0,kFALSE);
    TComplex cDenomNeg = TComplex(0.0,0.0,kFALSE);

    cDenom = Four(0,0,0,0);
    cNom = Four(fiHarm[0],fiHarm[1], fiHarm[2],fiHarm[3]);

    Bool_t bFillPos = kFALSE;
    Bool_t bFillNeg = kFALSE;
    
    
    if(cDenom.Re() > 0.0) { bFillPos = kTRUE; value = cNom.Re() / cDenom.Re(); }
    if(bFillPos && TMath::Abs(value) > 1.0) { bFillPos = kFALSE; }

    if(value > -100 && fCentrality< 90){
        Int_t random_number = rand() % 10 ;
        if(bFillPos) filling_collect[10]->Fill(fCentrality, value, cDenom.Re());
        if(bFillPos) filling_collect[random_number]->Fill(fCentrality, value, cDenom.Re());
    }
    
    return ;
    
}

TComplex AliAnalysisFlowTask5::Two(const Int_t n1, const Int_t n2) const
{
 TComplex formula = Q(n1,1)*Q(n2,1) - Q(n1+n2,2);
 return formula;
}

// ============================================================================
TComplex AliAnalysisFlowTask5::TwoGap(const Int_t n1, const Int_t n2) const
{
 TComplex formula = QGapPos(n1,1)*QGapNeg(n2,1);
 return formula;
}


TComplex AliAnalysisFlowTask5::Four(const Int_t n1, const Int_t n2, const Int_t n3, const Int_t n4) const
{
  TComplex formula = Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)-Q(n1+n2,2)*Q(n3,1)*Q(n4,1)-Q(n2,1)*Q(n1+n3,2)*Q(n4,1)
                    - Q(n1,1)*Q(n2+n3,2)*Q(n4,1)+2.0*Q(n1+n2+n3,3)*Q(n4,1)-Q(n2,1)*Q(n3,1)*Q(n1+n4,2)
                    + Q(n2+n3,2)*Q(n1+n4,2)-Q(n1,1)*Q(n3,1)*Q(n2+n4,2)+Q(n1+n3,2)*Q(n2+n4,2)
                    + 2.0*Q(n3,1)*Q(n1+n2+n4,3)-Q(n1,1)*Q(n2,1)*Q(n3+n4,2)+Q(n1+n2,2)*Q(n3+n4,2)
                    + 2.0*Q(n2,1)*Q(n1+n3+n4,3)+2.0*Q(n1,1)*Q(n2+n3+n4,3)-6.0*Q(n1+n2+n3+n4,4);
  return formula;
}


TComplex AliAnalysisFlowTask5::Three(const Int_t n1, const Int_t n2, const Int_t n3) const
{
  TComplex formula = Q(n1,1)*Q(n2,1)*Q(n3,1)-Q(n1+n2,2)*Q(n3,1)-Q(n2,1)*Q(n1+n3,2)
                          - Q(n1,1)*Q(n2+n3,2)+2.*Q(n1+n2+n3,3);
  return formula;
}
