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
    
    const char *name_SC_5_2[]= {"4pc_5_2_1", "4pc_5_2_2", "4pc_5_2_3", "4pc_5_2_4","4pc_5_2_5", "4pc_5_2_6", "4pc_5_2_7", "4pc_5_2_8", "4pc_5_2_9" ,"4pc_5_2_10"}; // enough to hold all numbers up to 64-bits TODO: FIX THIS FUCKIN SHIT
    const char *name_SC_5_3[]= {"4pc_5_3_1", "4pc_5_3_2", "4pc_5_3_3", "4pc_5_3_4","4pc_5_3_5", "4pc_5_3_6", "4pc_5_3_7", "4pc_5_3_8", "4pc_5_3_9" ,"4pc_5_3_10"}; // enough to hold all numbers up to 64-bits TODO: FIX THIS FUCKIN SHIT
    const char *name_SC_4_3[]= {"4pc_4_3_1", "4pc_4_3_2", "4pc_4_3_3", "4pc_4_3_4","4pc_4_3_5", "4pc_4_3_6", "4pc_4_3_7", "4pc_4_3_8", "4pc_4_3_9" ,"4pc_4_3_10"}; // enough to hold all numbers up to 64-bits TODO: FIX  FUCKIN SHIT
    const char *name_SC_4_2[]= {"4pc_4_2_1", "4pc_4_2_2", "4pc_4_2_3", "4pc_4_2_4","4pc_4_2_5", "4pc_4_2_6", "4pc_4_2_7", "4pc_4_2_8", "4pc_4_2_9" ,"4pc_4_2_10"}; // enough to hold all numbers up to 64-bits TODO: FIX  FUCKIN SHIT
    const char *name_SC_3_2[]= {"4pc_3_2_1", "4pc_3_2_2", "4pc_3_2_3", "4pc_3_2_4","4pc_3_2_5", "4pc_3_2_6", "4pc_3_2_7", "4pc_3_2_8", "4pc_3_2_9" ,"4pc_3_2_10"}; // enough to hold all numbers up to 64-bits TODO: FIX  FUCKIN SHIT

    const char *name_2pc_2[]= {"2pc_2_1", "2pc_2_2", "2pc_2_3", "2pc_2_4","2pc_2_5", "2pc_2_6", "2pc_2_7", "2pc_2_8", "2pc_2_9" ,"2pc_2_10"}; // enough to hold all numbers up to 64-bits TODO: FIX  FUCKIN SHIT
    const char *name_2pc_3[]= {"2pc_3_1", "2pc_3_2", "2pc_3_3", "2pc_3_4","2pc_3_5", "2pc_3_6", "2pc_3_7", "2pc_3_8", "2pc_3_9" ,"2pc_3_10"}; // enough to hold all numbers up to 64-bits TODO: FIX  FUCKIN SHIT
    const char *name_2pc_4[]= {"2pc_4_1", "2pc_4_2", "2pc_4_3", "2pc_4_4","2pc_4_5", "2pc_4_6", "2pc_4_7", "2pc_4_8", "2pc_4_9" ,"2pc_4_10"}; // enough to hold all numbers up to 64-bits TODO: FIX  FUCKIN SHIT
    const char *name_2pc_5[]= {"2pc_5_1", "2pc_5_2", "2pc_5_3", "2pc_5_4","2pc_5_5", "2pc_5_6", "2pc_5_7", "2pc_5_8", "2pc_5_9" ,"2pc_5_10"}; // enough to hold all numbers up to 64-bits TODO: FIX  FUCKIN SHIT
    
    fProfile_SC_5_2 = new TProfile("4_pc_5_2","Profile eta vs value; centrality;  <4>(5,2) ",6, bin_edges);
    fProfile_SC_5_3 = new TProfile("4_pc_5_3","Profile eta vs value; centrality;  <4>(5,3) ",6, bin_edges);
    fProfile_SC_4_3 = new TProfile("4_pc_4_3 ","Profile eta vs value; centrality;  <4>(4,3) ",6, bin_edges);
    fProfile_SC_4_2 = new TProfile("4_pc_4_2 ","Profile eta vs value; centrality;  <4>(4,2) ",6, bin_edges);
    fProfile_SC_3_2 = new TProfile("4_pc_3_2 ","Profile eta vs value; centrality;  <4>(3,2) ",6, bin_edges);

    fProfile_2pc_2 = new TProfile("2_pc_2 ","Profile eta vs value; centrality;  <2>(2) ",6, bin_edges);
    fProfile_2pc_3 = new TProfile("2_pc_3 ","Profile eta vs value; centrality;  <2>(3) ",6, bin_edges);
    fProfile_2pc_4 = new TProfile("2_pc_4 ","Profile eta vs value; centrality;  <2>(4) ",6, bin_edges);
    fProfile_2pc_5 = new TProfile("2_pc_5 ","Profile eta vs value; centrality;  <2>(5) ",6, bin_edges);
    
    fOutputList->Add(fProfile_SC_5_2);
    fOutputList->Add(fProfile_SC_5_3);
    fOutputList->Add(fProfile_SC_4_3);
    fOutputList->Add(fProfile_SC_4_2);
    fOutputList->Add(fProfile_SC_3_2);
    
    fOutputList->Add(fProfile_2pc_2);
    fOutputList->Add(fProfile_2pc_3);
    fOutputList->Add(fProfile_2pc_4);
    fOutputList->Add(fProfile_2pc_5);
    
    for(Int_t i(0); i<10; i++){
        fProfiles_collection_SC_5_2[i]= new TProfile(name_SC_5_2[i],"Profile eta vs value; centrality;  <4>(5,2) ",6, bin_edges);
        fProfiles_collection_SC_5_3[i]= new TProfile(name_SC_5_3[i],"Profile eta vs value; centrality;  <4>(5,3) ",6, bin_edges);
        fProfiles_collection_SC_4_3[i]= new TProfile(name_SC_4_3[i],"Profile eta vs value; centrality;  <4>(4,3) ",6, bin_edges);
        fProfiles_collection_SC_4_2[i]= new TProfile(name_SC_4_2[i],"Profile eta vs value; centrality;  <4>(4,2) ",6, bin_edges);
        fProfiles_collection_SC_3_2[i]= new TProfile(name_SC_3_2[i],"Profile eta vs value; centrality;  <4>(3,2) ",6, bin_edges);
        
        fProfiles_collection_2pc_2[i] = new TProfile(name_2pc_2[i],"Profile eta vs value; centrality;  <2>(2) ",6, bin_edges);
        fProfiles_collection_2pc_3[i] = new TProfile(name_2pc_3[i],"Profile eta vs value; centrality;  <2>(2) ",6, bin_edges);
        fProfiles_collection_2pc_4[i] = new TProfile(name_2pc_4[i],"Profile eta vs value; centrality;  <2>(2) ",6, bin_edges);
        fProfiles_collection_2pc_5[i] = new TProfile(name_2pc_5[i],"Profile eta vs value; centrality;  <2>(2) ",6, bin_edges);
        
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
    
    Calculate_correlation_4(fProfile_SC_5_2, fProfiles_collection_SC_5_2);
    
    fiHarm[0]=3;fiHarm[1]=5;fiHarm[2]=-5;fiHarm[3]=-3;
    
    Calculate_correlation_4(fProfile_SC_5_3, fProfiles_collection_SC_5_3);
    
    fiHarm[0]=3;fiHarm[1]=4;fiHarm[2]=-4;fiHarm[3]=-3;
    
    Calculate_correlation_4(fProfile_SC_4_3, fProfiles_collection_SC_4_3);
    
    fiHarm[0]=2;fiHarm[1]=4;fiHarm[2]=-4;fiHarm[3]=-2;
    
    Calculate_correlation_4(fProfile_SC_4_2, fProfiles_collection_SC_4_2);
    
    fiHarm[0]=2;fiHarm[1]=3;fiHarm[2]=-2;fiHarm[3]=-3;
    
    Calculate_correlation_4(fProfile_SC_3_2, fProfiles_collection_SC_3_2);
    
    fiHarm[0]=2;fiHarm[1]=-2;
    
    Calculate_correlation_2(fProfile_2pc_2, fProfiles_collection_2pc_2);
    
    fiHarm[0]=3;fiHarm[1]=-3;
    
    Calculate_correlation_2(fProfile_2pc_3, fProfiles_collection_2pc_3);
    
    fiHarm[0]=4;fiHarm[1]=-4;
    
    Calculate_correlation_2(fProfile_2pc_4, fProfiles_collection_2pc_4);
    
    fiHarm[0]=5;fiHarm[1]=-5;
    
    Calculate_correlation_2(fProfile_2pc_5, fProfiles_collection_2pc_5);
    
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
    /*
    if(track->DCA() > fDCAtoVtxXY) { return kFALSE; }
    if(track->ZAtDCA() > fDCAtoVtxZ) { return kFALSE; }
    if(track->GetTPCNcls() < fNTPCClusters)  { return kFALSE; }
    if(track->GetITSNcls()< fNITSClusters)  { return kFALSE; }
    if(track->GetTPCchi2perCluster() < fTPCchi2perClusterMin || track->GetTPCchi2perCluster() > fTPCchi2perClusterMax)  { return kFALSE; }
    */
    
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
/*
TComplex AliAnalysisFlowTask5::QGapNeg(const Int_t n, const Int_t p) const
{
 if(n < 0) return TComplex::Conjugate(fFlowVecQneg[-n][p]);
 else return fFlowVecQneg[n][p];
}
*/
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


void AliAnalysisFlowTask5::Calculate_correlation_2(TProfile* filling_prof, TProfile*  filling_collect[]){
    

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
        if(bFillPos) filling_prof ->Fill(fCentrality, value, cDenom.Re());
        if(bFillPos) filling_collect[random_number]->Fill(fCentrality, value, cDenom.Re());
        //!fMyProfile_v_SC_5_2_2->Fill(fCentrality,value, cDenom.Re());
    }
    
    return ;
    
}


void AliAnalysisFlowTask5::Calculate_correlation_4(TProfile* filling_prof, TProfile*  filling_collect[]){
    

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
        if(bFillPos) filling_prof ->Fill(fCentrality, value, cDenom.Re());
        if(bFillPos) filling_collect[random_number]->Fill(fCentrality, value, cDenom.Re());
    }
    
    return ;
    
}

/*
void AliAnalysisFlowTask5::Calculate_correlation_6(){
    

   //! cout<<cNom<<"\t";
    // important to distinguish between different "setup"
    
    Double_t value=0.;

    TComplex cNom = TComplex(0.0,0.0,kFALSE);
    TComplex cDenom = TComplex(0.0,0.0,kFALSE);
    TComplex cNomNeg = TComplex(0.0,0.0,kFALSE);
    TComplex cDenomNeg = TComplex(0.0,0.0,kFALSE);
    
    if(fDiff) {
        cDenom = TwoDiff(0,0);
        cNom = TwoDiff(fiHarm[0],fiHarm[1]);
    }

    else {
        cDenom = Six(0,0,0,0,0,0);
        cNom = Six(fiHarm[0],fiHarm[1],fiHarm[2],fiHarm[3],fiHarm[4],fiHarm[5]);
    }
    
    Bool_t bFillPos = kFALSE;
    Bool_t bFillNeg = kFALSE;
    
    
    if(cDenom.Re() > 0.0) { bFillPos = kTRUE; value = cNom.Re() / cDenom.Re(); }
    if(bFillPos && TMath::Abs(value) > 1.0) { bFillPos = kFALSE; }

    if(value > -10 && fCentrality< 90){
        Int_t random_number = rand() % 10 ;
        if(bFillPos) fMyProfile_v_SC_5_2_6->Fill(fCentrality, value, cDenom.Re());
        if(bFillPos) fProfiles_collection6[random_number]->Fill(fCentrality, value, cDenom.Re());
    }
    
    return ;
    
}

 */
/*
void AliAnalysisFlowTask5::Calculate_correlation_8(){
    

   //! cout<<cNom<<"\t";
    // important to distinguish between different "setup"
    
    Double_t value=0.;

    TComplex cNom = TComplex(0.0,0.0,kFALSE);
    TComplex cDenom = TComplex(0.0,0.0,kFALSE);
    TComplex cNomNeg = TComplex(0.0,0.0,kFALSE);
    TComplex cDenomNeg = TComplex(0.0,0.0,kFALSE);
    
    cDenom = Eight(0,0,0,0,0,0,0,0);
    cNom = Eight(fiHarm[0],fiHarm[1],fiHarm[2],fiHarm[3],fiHarm[4],fiHarm[5], fiHarm[6], fiHarm[7]);
    
    Bool_t bFillPos = kFALSE;
    Bool_t bFillNeg = kFALSE;
    
    if(cDenom.Re() > 0.0) { bFillPos = kTRUE; value = cNom.Re() / cDenom.Re(); }
    if(bFillPos && TMath::Abs(value) > 1.0) { bFillPos = kFALSE; }

    if(value > -10 && fCentrality< 90){
        Int_t random_number = rand() % 10 ;
        if(bFillPos) fMyProfile_v_SC_5_2_8->Fill(fCentrality, value, cDenom.Re());
        if(bFillPos) fProfiles_collection8[random_number]->Fill(fCentrality, value, cDenom.Re());
    }
    
    return ;
    
}


*/

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


TComplex AliAnalysisFlowTask5::Six(const Int_t n1, const Int_t n2, const Int_t n3, const Int_t n4, const Int_t n5, const Int_t n6) const
{
  TComplex formula = Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n6,1)-Q(n1+n2,2)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n6,1)
    - Q(n2,1)*Q(n1+n3,2)*Q(n4,1)*Q(n5,1)*Q(n6,1)-Q(n1,1)*Q(n2+n3,2)*Q(n4,1)*Q(n5,1)*Q(n6,1)
    + 2.*Q(n1+n2+n3,3)*Q(n4,1)*Q(n5,1)*Q(n6,1)-Q(n2,1)*Q(n3,1)*Q(n1+n4,2)*Q(n5,1)*Q(n6,1)
    + Q(n2+n3,2)*Q(n1+n4,2)*Q(n5,1)*Q(n6,1)-Q(n1,1)*Q(n3,1)*Q(n2+n4,2)*Q(n5,1)*Q(n6,1)
    + Q(n1+n3,2)*Q(n2+n4,2)*Q(n5,1)*Q(n6,1)+2.*Q(n3,1)*Q(n1+n2+n4,3)*Q(n5,1)*Q(n6,1)
    - Q(n1,1)*Q(n2,1)*Q(n3+n4,2)*Q(n5,1)*Q(n6,1)+Q(n1+n2,2)*Q(n3+n4,2)*Q(n5,1)*Q(n6,1)
    + 2.*Q(n2,1)*Q(n1+n3+n4,3)*Q(n5,1)*Q(n6,1)+2.*Q(n1,1)*Q(n2+n3+n4,3)*Q(n5,1)*Q(n6,1)
    - 6.*Q(n1+n2+n3+n4,4)*Q(n5,1)*Q(n6,1)-Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n1+n5,2)*Q(n6,1)
    + Q(n2+n3,2)*Q(n4,1)*Q(n1+n5,2)*Q(n6,1)+Q(n3,1)*Q(n2+n4,2)*Q(n1+n5,2)*Q(n6,1)
    + Q(n2,1)*Q(n3+n4,2)*Q(n1+n5,2)*Q(n6,1)-2.*Q(n2+n3+n4,3)*Q(n1+n5,2)*Q(n6,1)
    - Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n2+n5,2)*Q(n6,1)+Q(n1+n3,2)*Q(n4,1)*Q(n2+n5,2)*Q(n6,1)
    + Q(n3,1)*Q(n1+n4,2)*Q(n2+n5,2)*Q(n6,1)+Q(n1,1)*Q(n3+n4,2)*Q(n2+n5,2)*Q(n6,1)
    - 2.*Q(n1+n3+n4,3)*Q(n2+n5,2)*Q(n6,1)+2.*Q(n3,1)*Q(n4,1)*Q(n1+n2+n5,3)*Q(n6,1)
    - 2.*Q(n3+n4,2)*Q(n1+n2+n5,3)*Q(n6,1)-Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n3+n5,2)*Q(n6,1)
    + Q(n1+n2,2)*Q(n4,1)*Q(n3+n5,2)*Q(n6,1)+Q(n2,1)*Q(n1+n4,2)*Q(n3+n5,2)*Q(n6,1)
    + Q(n1,1)*Q(n2+n4,2)*Q(n3+n5,2)*Q(n6,1)-2.*Q(n1+n2+n4,3)*Q(n3+n5,2)*Q(n6,1)
    + 2.*Q(n2,1)*Q(n4,1)*Q(n1+n3+n5,3)*Q(n6,1)-2.*Q(n2+n4,2)*Q(n1+n3+n5,3)*Q(n6,1)
    + 2.*Q(n1,1)*Q(n4,1)*Q(n2+n3+n5,3)*Q(n6,1)-2.*Q(n1+n4,2)*Q(n2+n3+n5,3)*Q(n6,1)
    - 6.*Q(n4,1)*Q(n1+n2+n3+n5,4)*Q(n6,1)-Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4+n5,2)*Q(n6,1)
    + Q(n1+n2,2)*Q(n3,1)*Q(n4+n5,2)*Q(n6,1)+Q(n2,1)*Q(n1+n3,2)*Q(n4+n5,2)*Q(n6,1)
    + Q(n1,1)*Q(n2+n3,2)*Q(n4+n5,2)*Q(n6,1)-2.*Q(n1+n2+n3,3)*Q(n4+n5,2)*Q(n6,1)
    + 2.*Q(n2,1)*Q(n3,1)*Q(n1+n4+n5,3)*Q(n6,1)-2.*Q(n2+n3,2)*Q(n1+n4+n5,3)*Q(n6,1)
    + 2.*Q(n1,1)*Q(n3,1)*Q(n2+n4+n5,3)*Q(n6,1)-2.*Q(n1+n3,2)*Q(n2+n4+n5,3)*Q(n6,1)
    - 6.*Q(n3,1)*Q(n1+n2+n4+n5,4)*Q(n6,1)+2.*Q(n1,1)*Q(n2,1)*Q(n3+n4+n5,3)*Q(n6,1)
    - 2.*Q(n1+n2,2)*Q(n3+n4+n5,3)*Q(n6,1)-6.*Q(n2,1)*Q(n1+n3+n4+n5,4)*Q(n6,1)
    - 6.*Q(n1,1)*Q(n2+n3+n4+n5,4)*Q(n6,1)+24.*Q(n1+n2+n3+n4+n5,5)*Q(n6,1)
    - Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n1+n6,2)+Q(n2+n3,2)*Q(n4,1)*Q(n5,1)*Q(n1+n6,2)
    + Q(n3,1)*Q(n2+n4,2)*Q(n5,1)*Q(n1+n6,2)+Q(n2,1)*Q(n3+n4,2)*Q(n5,1)*Q(n1+n6,2)
    - 2.*Q(n2+n3+n4,3)*Q(n5,1)*Q(n1+n6,2)+Q(n3,1)*Q(n4,1)*Q(n2+n5,2)*Q(n1+n6,2)
    - Q(n3+n4,2)*Q(n2+n5,2)*Q(n1+n6,2)+Q(n2,1)*Q(n4,1)*Q(n3+n5,2)*Q(n1+n6,2)
    - Q(n2+n4,2)*Q(n3+n5,2)*Q(n1+n6,2)-2.*Q(n4,1)*Q(n2+n3+n5,3)*Q(n1+n6,2)
    + Q(n2,1)*Q(n3,1)*Q(n4+n5,2)*Q(n1+n6,2)-Q(n2+n3,2)*Q(n4+n5,2)*Q(n1+n6,2)
    - 2.*Q(n3,1)*Q(n2+n4+n5,3)*Q(n1+n6,2)-2.*Q(n2,1)*Q(n3+n4+n5,3)*Q(n1+n6,2)
    + 6.*Q(n2+n3+n4+n5,4)*Q(n1+n6,2)-Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n2+n6,2)
    + Q(n1+n3,2)*Q(n4,1)*Q(n5,1)*Q(n2+n6,2)+Q(n3,1)*Q(n1+n4,2)*Q(n5,1)*Q(n2+n6,2)
    + Q(n1,1)*Q(n3+n4,2)*Q(n5,1)*Q(n2+n6,2)-2.*Q(n1+n3+n4,3)*Q(n5,1)*Q(n2+n6,2)
    + Q(n3,1)*Q(n4,1)*Q(n1+n5,2)*Q(n2+n6,2)-Q(n3+n4,2)*Q(n1+n5,2)*Q(n2+n6,2)
    + Q(n1,1)*Q(n4,1)*Q(n3+n5,2)*Q(n2+n6,2)-Q(n1+n4,2)*Q(n3+n5,2)*Q(n2+n6,2)
    - 2.*Q(n4,1)*Q(n1+n3+n5,3)*Q(n2+n6,2)+Q(n1,1)*Q(n3,1)*Q(n4+n5,2)*Q(n2+n6,2)
    - Q(n1+n3,2)*Q(n4+n5,2)*Q(n2+n6,2)-2.*Q(n3,1)*Q(n1+n4+n5,3)*Q(n2+n6,2)
    - 2.*Q(n1,1)*Q(n3+n4+n5,3)*Q(n2+n6,2)+6.*Q(n1+n3+n4+n5,4)*Q(n2+n6,2)
    + 2.*Q(n3,1)*Q(n4,1)*Q(n5,1)*Q(n1+n2+n6,3)-2.*Q(n3+n4,2)*Q(n5,1)*Q(n1+n2+n6,3)
    - 2.*Q(n4,1)*Q(n3+n5,2)*Q(n1+n2+n6,3)-2.*Q(n3,1)*Q(n4+n5,2)*Q(n1+n2+n6,3)
    + 4.*Q(n3+n4+n5,3)*Q(n1+n2+n6,3)-Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n5,1)*Q(n3+n6,2)
    + Q(n1+n2,2)*Q(n4,1)*Q(n5,1)*Q(n3+n6,2)+Q(n2,1)*Q(n1+n4,2)*Q(n5,1)*Q(n3+n6,2)
    + Q(n1,1)*Q(n2+n4,2)*Q(n5,1)*Q(n3+n6,2)-2.*Q(n1+n2+n4,3)*Q(n5,1)*Q(n3+n6,2)
    + Q(n2,1)*Q(n4,1)*Q(n1+n5,2)*Q(n3+n6,2)-Q(n2+n4,2)*Q(n1+n5,2)*Q(n3+n6,2)
    + Q(n1,1)*Q(n4,1)*Q(n2+n5,2)*Q(n3+n6,2)-Q(n1+n4,2)*Q(n2+n5,2)*Q(n3+n6,2)
    - 2.*Q(n4,1)*Q(n1+n2+n5,3)*Q(n3+n6,2)+Q(n1,1)*Q(n2,1)*Q(n4+n5,2)*Q(n3+n6,2)
    - Q(n1+n2,2)*Q(n4+n5,2)*Q(n3+n6,2)-2.*Q(n2,1)*Q(n1+n4+n5,3)*Q(n3+n6,2)
    - 2.*Q(n1,1)*Q(n2+n4+n5,3)*Q(n3+n6,2)+6.*Q(n1+n2+n4+n5,4)*Q(n3+n6,2)
    + 2.*Q(n2,1)*Q(n4,1)*Q(n5,1)*Q(n1+n3+n6,3)-2.*Q(n2+n4,2)*Q(n5,1)*Q(n1+n3+n6,3)
    - 2.*Q(n4,1)*Q(n2+n5,2)*Q(n1+n3+n6,3)-2.*Q(n2,1)*Q(n4+n5,2)*Q(n1+n3+n6,3)
    + 4.*Q(n2+n4+n5,3)*Q(n1+n3+n6,3)+2.*Q(n1,1)*Q(n4,1)*Q(n5,1)*Q(n2+n3+n6,3)
    - 2.*Q(n1+n4,2)*Q(n5,1)*Q(n2+n3+n6,3)-2.*Q(n4,1)*Q(n1+n5,2)*Q(n2+n3+n6,3)
    - 2.*Q(n1,1)*Q(n4+n5,2)*Q(n2+n3+n6,3)+4.*Q(n1+n4+n5,3)*Q(n2+n3+n6,3)
    - 6.*Q(n4,1)*Q(n5,1)*Q(n1+n2+n3+n6,4)+6.*Q(n4+n5,2)*Q(n1+n2+n3+n6,4)
    - Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n5,1)*Q(n4+n6,2)+Q(n1+n2,2)*Q(n3,1)*Q(n5,1)*Q(n4+n6,2)
    + Q(n2,1)*Q(n1+n3,2)*Q(n5,1)*Q(n4+n6,2)+Q(n1,1)*Q(n2+n3,2)*Q(n5,1)*Q(n4+n6,2)
    - 2.*Q(n1+n2+n3,3)*Q(n5,1)*Q(n4+n6,2)+Q(n2,1)*Q(n3,1)*Q(n1+n5,2)*Q(n4+n6,2)
    - Q(n2+n3,2)*Q(n1+n5,2)*Q(n4+n6,2)+Q(n1,1)*Q(n3,1)*Q(n2+n5,2)*Q(n4+n6,2)
    - Q(n1+n3,2)*Q(n2+n5,2)*Q(n4+n6,2)-2.*Q(n3,1)*Q(n1+n2+n5,3)*Q(n4+n6,2)
    + Q(n1,1)*Q(n2,1)*Q(n3+n5,2)*Q(n4+n6,2)-Q(n1+n2,2)*Q(n3+n5,2)*Q(n4+n6,2)
    - 2.*Q(n2,1)*Q(n1+n3+n5,3)*Q(n4+n6,2)-2.*Q(n1,1)*Q(n2+n3+n5,3)*Q(n4+n6,2)
    + 6.*Q(n1+n2+n3+n5,4)*Q(n4+n6,2)+2.*Q(n2,1)*Q(n3,1)*Q(n5,1)*Q(n1+n4+n6,3)
    - 2.*Q(n2+n3,2)*Q(n5,1)*Q(n1+n4+n6,3)-2.*Q(n3,1)*Q(n2+n5,2)*Q(n1+n4+n6,3)
    - 2.*Q(n2,1)*Q(n3+n5,2)*Q(n1+n4+n6,3)+4.*Q(n2+n3+n5,3)*Q(n1+n4+n6,3)
    + 2.*Q(n1,1)*Q(n3,1)*Q(n5,1)*Q(n2+n4+n6,3)-2.*Q(n1+n3,2)*Q(n5,1)*Q(n2+n4+n6,3)
    - 2.*Q(n3,1)*Q(n1+n5,2)*Q(n2+n4+n6,3)-2.*Q(n1,1)*Q(n3+n5,2)*Q(n2+n4+n6,3)
    + 4.*Q(n1+n3+n5,3)*Q(n2+n4+n6,3)-6.*Q(n3,1)*Q(n5,1)*Q(n1+n2+n4+n6,4)
    + 6.*Q(n3+n5,2)*Q(n1+n2+n4+n6,4)+2.*Q(n1,1)*Q(n2,1)*Q(n5,1)*Q(n3+n4+n6,3)
    - 2.*Q(n1+n2,2)*Q(n5,1)*Q(n3+n4+n6,3)-2.*Q(n2,1)*Q(n1+n5,2)*Q(n3+n4+n6,3)
    - 2.*Q(n1,1)*Q(n2+n5,2)*Q(n3+n4+n6,3)+4.*Q(n1+n2+n5,3)*Q(n3+n4+n6,3)
    - 6.*Q(n2,1)*Q(n5,1)*Q(n1+n3+n4+n6,4)+6.*Q(n2+n5,2)*Q(n1+n3+n4+n6,4)
    - 6.*Q(n1,1)*Q(n5,1)*Q(n2+n3+n4+n6,4)+6.*Q(n1+n5,2)*Q(n2+n3+n4+n6,4)
    + 24.*Q(n5,1)*Q(n1+n2+n3+n4+n6,5)-Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5+n6,2)
    + Q(n1+n2,2)*Q(n3,1)*Q(n4,1)*Q(n5+n6,2)+Q(n2,1)*Q(n1+n3,2)*Q(n4,1)*Q(n5+n6,2)
    + Q(n1,1)*Q(n2+n3,2)*Q(n4,1)*Q(n5+n6,2)-2.*Q(n1+n2+n3,3)*Q(n4,1)*Q(n5+n6,2)
    + Q(n2,1)*Q(n3,1)*Q(n1+n4,2)*Q(n5+n6,2)-Q(n2+n3,2)*Q(n1+n4,2)*Q(n5+n6,2)
    + Q(n1,1)*Q(n3,1)*Q(n2+n4,2)*Q(n5+n6,2)-Q(n1+n3,2)*Q(n2+n4,2)*Q(n5+n6,2)
    - 2.*Q(n3,1)*Q(n1+n2+n4,3)*Q(n5+n6,2)+Q(n1,1)*Q(n2,1)*Q(n3+n4,2)*Q(n5+n6,2)
    - Q(n1+n2,2)*Q(n3+n4,2)*Q(n5+n6,2)-2.*Q(n2,1)*Q(n1+n3+n4,3)*Q(n5+n6,2)
    - 2.*Q(n1,1)*Q(n2+n3+n4,3)*Q(n5+n6,2)+6.*Q(n1+n2+n3+n4,4)*Q(n5+n6,2)
    + 2.*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n1+n5+n6,3)-2.*Q(n2+n3,2)*Q(n4,1)*Q(n1+n5+n6,3)
    - 2.*Q(n3,1)*Q(n2+n4,2)*Q(n1+n5+n6,3)-2.*Q(n2,1)*Q(n3+n4,2)*Q(n1+n5+n6,3)
    + 4.*Q(n2+n3+n4,3)*Q(n1+n5+n6,3)+2.*Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n2+n5+n6,3)
    - 2.*Q(n1+n3,2)*Q(n4,1)*Q(n2+n5+n6,3)-2.*Q(n3,1)*Q(n1+n4,2)*Q(n2+n5+n6,3)
    - 2.*Q(n1,1)*Q(n3+n4,2)*Q(n2+n5+n6,3)+4.*Q(n1+n3+n4,3)*Q(n2+n5+n6,3)
    - 6.*Q(n3,1)*Q(n4,1)*Q(n1+n2+n5+n6,4)+6.*Q(n3+n4,2)*Q(n1+n2+n5+n6,4)
    + 2.*Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n3+n5+n6,3)-2.*Q(n1+n2,2)*Q(n4,1)*Q(n3+n5+n6,3)
    - 2.*Q(n2,1)*Q(n1+n4,2)*Q(n3+n5+n6,3)-2.*Q(n1,1)*Q(n2+n4,2)*Q(n3+n5+n6,3)
    + 4.*Q(n1+n2+n4,3)*Q(n3+n5+n6,3)-6.*Q(n2,1)*Q(n4,1)*Q(n1+n3+n5+n6,4)
    + 6.*Q(n2+n4,2)*Q(n1+n3+n5+n6,4)-6.*Q(n1,1)*Q(n4,1)*Q(n2+n3+n5+n6,4)
    + 6.*Q(n1+n4,2)*Q(n2+n3+n5+n6,4)+24.*Q(n4,1)*Q(n1+n2+n3+n5+n6,5)
    + 2.*Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4+n5+n6,3)-2.*Q(n1+n2,2)*Q(n3,1)*Q(n4+n5+n6,3)
    - 2.*Q(n2,1)*Q(n1+n3,2)*Q(n4+n5+n6,3)-2.*Q(n1,1)*Q(n2+n3,2)*Q(n4+n5+n6,3)
    + 4.*Q(n1+n2+n3,3)*Q(n4+n5+n6,3)-6.*Q(n2,1)*Q(n3,1)*Q(n1+n4+n5+n6,4)
    + 6.*Q(n2+n3,2)*Q(n1+n4+n5+n6,4)-6.*Q(n1,1)*Q(n3,1)*Q(n2+n4+n5+n6,4)
    + 6.*Q(n1+n3,2)*Q(n2+n4+n5+n6,4)+24.*Q(n3,1)*Q(n1+n2+n4+n5+n6,5)
    - 6.*Q(n1,1)*Q(n2,1)*Q(n3+n4+n5+n6,4)+6.*Q(n1+n2,2)*Q(n3+n4+n5+n6,4)
    + 24.*Q(n2,1)*Q(n1+n3+n4+n5+n6,5)+24.*Q(n1,1)*Q(n2+n3+n4+n5+n6,5)
    - 120.*Q(n1+n2+n3+n4+n5+n6,6);
    return formula;
}

TComplex AliAnalysisFlowTask5::Eight(const Int_t n1, const Int_t n2, const Int_t n3, const Int_t n4, const Int_t n5, const Int_t n6, const Int_t n7, const Int_t n8) const
{
  TComplex Correlation = {0, 0};
    Int_t Narray[] = {n1, n2, n3, n4, n5, n6, n7};

    for(Int_t k=8; k-->0; )
    {// backward loop of k from m-1 until 0, where m is the m-particle correlation, in this case m=4

        Int_t array[7] = {0,1,2,3,4,5,6};
        Int_t iPerm = 0;
        Int_t count = 0;

        // k==7: there is just one combination, we can add it manually
        if(k==7){
            Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*
            Seven(n1, n2, n3, n4, n5, n6, n7)*Q(n8, 8-k);
        }// k==7

        else if(k==6){
            do{
                iPerm += 1;
                if(array[0] < array[1] && array[1] < array[2] && array[2] < array[3] && array[3] < array[4] && array[4] < array[5]){
                    count += 1;
                    Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*
                    Six(Narray[Int_t(array[0])], Narray[Int_t(array[1])], Narray[Int_t(array[2])],
                        Narray[Int_t(array[3])], Narray[Int_t(array[4])], Narray[Int_t(array[5])])*
                    Q(Narray[Int_t(array[6])]+n8, 8-k);
                }
            }while(std::next_permutation(array, array+7));
        }// k==6

        else if(k==5){
            do{
                iPerm += 1;
                if(iPerm%2 == 1){
                    if(array[0] < array[1] && array[1] < array[2] && array[2] < array[3] && array[3] < array[4]){
                        Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*
                        Five(Narray[Int_t(array[0])], Narray[Int_t(array[1])], Narray[Int_t(array[2])],
                             Narray[Int_t(array[3])], Narray[Int_t(array[4])])*
                        Q(Narray[Int_t(array[5])]+Narray[Int_t(array[6])]+n8, 8-k);
                    }
                }
            }while(std::next_permutation(array, array+7));
        }// k==5

        else if(k==4){
            do{
                iPerm += 1;
                if(iPerm%6 == 1){
                    if(array[0] < array[1] && array[1] < array[2] && array[2] < array[3]){
                        Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*
                        Four(Narray[Int_t(array[0])], Narray[Int_t(array[1])], Narray[Int_t(array[2])], Narray[Int_t(array[3])])*
                        Q(Narray[Int_t(array[4])]+Narray[Int_t(array[5])]+Narray[Int_t(array[6])]+n8, 8-k);
                    }
                }
            }while(std::next_permutation(array, array+7));
        }// k==4

        else if(k==3){
            do{
                iPerm += 1;
                if(iPerm%24 == 1){
                    if(array[0] < array[1] && array[1] < array[2]){
                        Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*
                        Three(Narray[Int_t(array[0])], Narray[Int_t(array[1])], Narray[Int_t(array[2])])*
                        Q(Narray[Int_t(array[3])]+Narray[Int_t(array[4])]+Narray[Int_t(array[5])]+Narray[Int_t(array[6])]+n8, 8-k);
                    }
                }
            }while(std::next_permutation(array, array+7));
        }// k==3

        else if(k==2){
            do{
                iPerm += 1;
                if(iPerm%120 == 1){
                    if(array[0] < array[1]){
                        Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*
                        Two(Narray[Int_t(array[0])], Narray[Int_t(array[1])])*
                        Q(Narray[Int_t(array[2])]+Narray[Int_t(array[3])]+Narray[Int_t(array[4])]
                          +Narray[Int_t(array[5])]+Narray[Int_t(array[6])]+n8, 8-k);
                    }
                }
            }while(std::next_permutation(array, array+7));
        }// k==2

        else if(k == 1){
            Correlation = Correlation
            + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*Q(n1, 1)*Q(n2+n3+n4+n5+n6+n7+n8, 8-k)
            + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*Q(n2, 1)*Q(n1+n3+n4+n5+n6+n7+n8, 8-k)
            + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*Q(n3, 1)*Q(n1+n2+n4+n5+n6+n7+n8, 8-k)
            + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*Q(n4, 1)*Q(n1+n2+n3+n5+n6+n7+n8, 8-k)
            + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*Q(n5, 1)*Q(n1+n2+n3+n4+n6+n7+n8, 8-k)
            + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*Q(n6, 1)*Q(n1+n2+n3+n4+n5+n7+n8, 8-k)
            + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*Q(n7, 1)*Q(n1+n2+n3+n4+n5+n6+n8, 8-k);
        }// k==1

        else if(k == 0){
            Correlation = Correlation + TMath::Power(-1, 8-k-1)*TMath::Factorial(8-k-1)*Q(n1+n2+n3+n4+n5+n6+n7+n8, 8-k);
        }// k==0

        else{
            printf("Invalid range of k in Eight() \n");
            return {0,0};
        }

    }// loop over k

    return Correlation;
}

TComplex AliAnalysisFlowTask5::Seven(const Int_t n1, const Int_t n2, const Int_t n3, const Int_t n4, const Int_t n5, const Int_t n6, const Int_t n7) const
{
  TComplex Correlation = {0, 0};
    Int_t Narray[] = {n1, n2, n3, n4, n5, n6};

    for(Int_t k=7; k-- >0; )
    {// backward loop of k from m-1 until 0, where m is the m-particle correlation, in this case m=4

        Int_t array[6] = {0,1,2,3,4,5};
        Int_t iPerm = 0;
        Int_t count = 0;

        // k==6: there is just one combination, we can add it manually
        if(k==6){
            Correlation = Correlation + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*
            Six(n1, n2, n3, n4, n5, n6)*Q(n7, 7-k);
        }// k==6

        else if(k==5){
            do{
                iPerm += 1;
                if(array[0] < array[1] && array[1] < array[2] && array[2] < array[3] && array[3] < array[4]){
                    count += 1;
                    Correlation = Correlation + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*
                    Five(Narray[Int_t(array[0])], Narray[Int_t(array[1])], Narray[Int_t(array[2])],
                         Narray[Int_t(array[3])], Narray[Int_t(array[4])])*
                    Q(Narray[int(array[5])]+n7, 7-k);
                }
            }while(std::next_permutation(array, array+6));
        }// k==5

        else if(k==4){
            do{
                iPerm += 1;
                if(iPerm%2 == 1){
                    if(array[0] < array[1] && array[1] < array[2] && array[2] < array[3]){
                        Correlation = Correlation + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*
                        Four(Narray[Int_t(array[0])], Narray[Int_t(array[1])], Narray[Int_t(array[2])],
                             Narray[Int_t(array[3])])*
                        Q(Narray[Int_t(array[4])]+Narray[Int_t(array[5])]+n7, 7-k);
                    }
                }
            }while(std::next_permutation(array, array+6));
        }// k==4

        else if(k==3){
            do{
                iPerm += 1;
                if(iPerm%6 == 1){
                    if(array[0] < array[1] && array[1] < array[2]){
                        Correlation = Correlation + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*
                        Three(Narray[Int_t(array[0])], Narray[Int_t(array[1])], Narray[Int_t(array[2])])*
                        Q(Narray[Int_t(array[3])]+Narray[Int_t(array[4])]+Narray[Int_t(array[5])]+n7, 7-k);
                    }
                }
            }while(std::next_permutation(array, array+6));
        }// k==3

        else if(k==2){
            do{
                iPerm += 1;
                if(iPerm%24 == 1){
                    if(array[0] < array[1]){
                        Correlation = Correlation + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*
                        Two(Narray[Int_t(array[0])], Narray[Int_t(array[1])])*
                        Q(Narray[Int_t(array[2])]+Narray[Int_t(array[3])]+Narray[Int_t(array[4])]
                          +Narray[Int_t(array[5])]+n7, 7-k);
                    }
                }
            }while(std::next_permutation(array, array+6));
        }// k==2

        else if(k == 1){
            Correlation = Correlation
            + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*Q(n1, 1)*Q(n2+n3+n4+n5+n6+n7, 7-k)
            + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*Q(n2, 1)*Q(n1+n3+n4+n5+n6+n7, 7-k)
            + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*Q(n3, 1)*Q(n1+n2+n4+n5+n6+n7, 7-k)
            + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*Q(n4, 1)*Q(n1+n2+n3+n5+n6+n7, 7-k)
            + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*Q(n5, 1)*Q(n1+n2+n3+n4+n6+n7, 7-k)
            + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*Q(n6, 1)*Q(n1+n2+n3+n4+n5+n7, 7-k);
        }// k==1

        else if(k == 0){
            Correlation = Correlation + TMath::Power(-1, 7-k-1)*TMath::Factorial(7-k-1)*Q(n1+n2+n3+n4+n5+n6+n7, 7-k);
        }// k==0

        else{
            printf("Invalid range of k in Seven()\n");;
            return {0,0};
        }

    }// loop over k
    return Correlation;
}
TComplex AliAnalysisFlowTask5::Five(const Int_t n1, const Int_t n2, const Int_t n3, const Int_t n4, const Int_t n5) const
{
    TComplex formula = Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n5,1)-Q(n1+n2,2)*Q(n3,1)*Q(n4,1)*Q(n5,1)
    - Q(n2,1)*Q(n1+n3,2)*Q(n4,1)*Q(n5,1)-Q(n1,1)*Q(n2+n3,2)*Q(n4,1)*Q(n5,1)
    + 2.*Q(n1+n2+n3,3)*Q(n4,1)*Q(n5,1)-Q(n2,1)*Q(n3,1)*Q(n1+n4,2)*Q(n5,1)
    + Q(n2+n3,2)*Q(n1+n4,2)*Q(n5,1)-Q(n1,1)*Q(n3,1)*Q(n2+n4,2)*Q(n5,1)
    + Q(n1+n3,2)*Q(n2+n4,2)*Q(n5,1)+2.*Q(n3,1)*Q(n1+n2+n4,3)*Q(n5,1)
    - Q(n1,1)*Q(n2,1)*Q(n3+n4,2)*Q(n5,1)+Q(n1+n2,2)*Q(n3+n4,2)*Q(n5,1)
    + 2.*Q(n2,1)*Q(n1+n3+n4,3)*Q(n5,1)+2.*Q(n1,1)*Q(n2+n3+n4,3)*Q(n5,1)
    - 6.*Q(n1+n2+n3+n4,4)*Q(n5,1)-Q(n2,1)*Q(n3,1)*Q(n4,1)*Q(n1+n5,2)
    + Q(n2+n3,2)*Q(n4,1)*Q(n1+n5,2)+Q(n3,1)*Q(n2+n4,2)*Q(n1+n5,2)
    + Q(n2,1)*Q(n3+n4,2)*Q(n1+n5,2)-2.*Q(n2+n3+n4,3)*Q(n1+n5,2)
    - Q(n1,1)*Q(n3,1)*Q(n4,1)*Q(n2+n5,2)+Q(n1+n3,2)*Q(n4,1)*Q(n2+n5,2)
    + Q(n3,1)*Q(n1+n4,2)*Q(n2+n5,2)+Q(n1,1)*Q(n3+n4,2)*Q(n2+n5,2)
    - 2.*Q(n1+n3+n4,3)*Q(n2+n5,2)+2.*Q(n3,1)*Q(n4,1)*Q(n1+n2+n5,3)
    - 2.*Q(n3+n4,2)*Q(n1+n2+n5,3)-Q(n1,1)*Q(n2,1)*Q(n4,1)*Q(n3+n5,2)
    + Q(n1+n2,2)*Q(n4,1)*Q(n3+n5,2)+Q(n2,1)*Q(n1+n4,2)*Q(n3+n5,2)
    + Q(n1,1)*Q(n2+n4,2)*Q(n3+n5,2)-2.*Q(n1+n2+n4,3)*Q(n3+n5,2)
    + 2.*Q(n2,1)*Q(n4,1)*Q(n1+n3+n5,3)-2.*Q(n2+n4,2)*Q(n1+n3+n5,3)
    + 2.*Q(n1,1)*Q(n4,1)*Q(n2+n3+n5,3)-2.*Q(n1+n4,2)*Q(n2+n3+n5,3)
    - 6.*Q(n4,1)*Q(n1+n2+n3+n5,4)-Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4+n5,2)
    + Q(n1+n2,2)*Q(n3,1)*Q(n4+n5,2)+Q(n2,1)*Q(n1+n3,2)*Q(n4+n5,2)
    + Q(n1,1)*Q(n2+n3,2)*Q(n4+n5,2)-2.*Q(n1+n2+n3,3)*Q(n4+n5,2)
    + 2.*Q(n2,1)*Q(n3,1)*Q(n1+n4+n5,3)-2.*Q(n2+n3,2)*Q(n1+n4+n5,3)
    + 2.*Q(n1,1)*Q(n3,1)*Q(n2+n4+n5,3)-2.*Q(n1+n3,2)*Q(n2+n4+n5,3)
    - 6.*Q(n3,1)*Q(n1+n2+n4+n5,4)+2.*Q(n1,1)*Q(n2,1)*Q(n3+n4+n5,3)
    - 2.*Q(n1+n2,2)*Q(n3+n4+n5,3)-6.*Q(n2,1)*Q(n1+n3+n4+n5,4)
    - 6.*Q(n1,1)*Q(n2+n3+n4+n5,4)+24.*Q(n1+n2+n3+n4+n5,5);
    return formula;
}

TComplex AliAnalysisFlowTask5::Three(const Int_t n1, const Int_t n2, const Int_t n3) const
{
  TComplex formula = Q(n1,1)*Q(n2,1)*Q(n3,1)-Q(n1+n2,2)*Q(n3,1)-Q(n2,1)*Q(n1+n3,2)
                          - Q(n1,1)*Q(n2+n3,2)+2.*Q(n1+n2+n3,3);
  return formula;
}
