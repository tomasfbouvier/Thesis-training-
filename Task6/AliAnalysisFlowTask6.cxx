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

/* AliAnalysisFlowTask6
 *
 * example flow class:
 * event selection
 * track selection
 * loading particle weights if needed (possible to change run (in)dependent)
 * filling Q vectors (only needed combinations of harmonics & powers! check if this is what you need!)
 * calculate m-particle correlations (using generic algorithm)
 */

#include "AliAnalysisFlowTask6.h"

class AliAnalysisFlowTask6;

using namespace std;

ClassImp(AliAnalysisFlowTask6)

AliAnalysisFlowTask6::AliAnalysisFlowTask6() : AliAnalysisTaskSE(),
    fAOD(0),
    fOutputListCharged(0),
    fTrigger(AliVEvent::kINT7),
    fIsPP(kFALSE),
    fEventRejectAddPileUp(kFALSE),
    fUseWeights(kFALSE),
    fFilterBit(96),
    fCentBinNum(10),
    fCentMin(0),
    fCentMax(100),
    fCentrality(0.0),
    fPtMin(0.2),
    fPtMax(5.0),
    fAbsEtaMin(0.4),
    fAbsEtaMax(0.8),
    fPVtxCutZ(10.0),
    fPtMinTrig(0.5),
    fPtMaxTrig(10.0),
    fPtMinAss(0.5),
    fPtMaxAss(1.5),
    fPVz(99.9),
    fCentEstimator("V0M"),
    fhEventCounter(0),
    fDCAtoVtxZ(3.0),
    fDCAtoVtxXY(3.0),
    fNTPCClusters(70),
    fNITSClusters(2),
    fTPCchi2perClusterMin(0.1),
    fTPCchi2perClusterMax(4.),
    fNzVtxBins(10),
    fNCentBins(15),
    fPoolMaxNEvents(2000),
    fPoolMinNTracks(50000),
    fEfficiencyEtaDependent(kFALSE),
    fUseEfficiency(kFALSE)
{}
//_____________________________________________________________________________
AliAnalysisFlowTask6::AliAnalysisFlowTask6(const char* name, Bool_t bUseEff) : AliAnalysisTaskSE(name),
    fAOD(0),
    fOutputListCharged(0),
    fTrigger(AliVEvent::kINT7),
    fIsPP(kFALSE),
    fEventRejectAddPileUp(kFALSE),
    fUseWeights(bUseEff),
    fFilterBit(96),
    fCentBinNum(10),
    fCentMin(0),
    fCentMax(100),
    fCentrality(0.0),
    fPtMin(0.2),
    fPtMax(5.0),
    fAbsEtaMin(0.4),
    fAbsEtaMax(0.8),
    fPVtxCutZ(10.0),
    fPtMinTrig(0.5),
    fPtMaxTrig(10.0),
    fPtMinAss(0.5),
    fPtMaxAss(1.5),
    fPVz(99.9),
    fCentEstimator("V0M"),
    fhEventCounter(0),
    fDCAtoVtxZ(3.0),
    fDCAtoVtxXY(3.0),
    fNTPCClusters(70),
    fNITSClusters(2),
    fTPCchi2perClusterMin(0.1),
    fTPCchi2perClusterMax(4.),
    fNzVtxBins(10),
    fNCentBins(15),
    fPoolMaxNEvents(2000),
    fPoolMinNTracks(50000),
    fEfficiencyEtaDependent(kFALSE),
    fUseEfficiency(bUseEff)

{
    DefineInput(0, TChain::Class());
    if(bUseEff) { DefineInput(1, TList::Class()); }
    DefineOutput(1, TList::Class());
}
//_____________________________________________________________________________
AliAnalysisFlowTask6::~AliAnalysisFlowTask6()
{
    if(fOutputListCharged) delete fOutputListCharged;

}
//_____________________________________________________________________________
void AliAnalysisFlowTask6::UserCreateOutputObjects()
{

    OpenFile(1);

    //!  ZUZANA'S ADD


    // //just for testing
    fPtBinsTrigCharged = {0.5, 1.0, 1.5, 2.0, 3.0, 5.0}; // Esto estaba muteado
    fPtBinsAss = {0.5, 1.0, 1.5, 2.0, 3.0};   //Esto estaba muteado
    fzVtxBins = {-10.0,-8.0,-6.0,-4.0,-2.0,0.0,2.0,4.0,6.0,8.0,10.0};
    fCentBins = {0,1,2,3,4,5,10,20,30,40,50,60,70,80,90,100};

    fOutputListCharged = new TList();
    fOutputListCharged->SetOwner(kTRUE);

    fhEventCounter = new TH1D("fhEventCounter","Event Counter",10,0,10);
    fOutputListCharged->Add(fhEventCounter);

    fHistPhiEta = new TH2D("fHistPhiEta", "fHistPhiEta; phi; eta", 100, -0.5, 7, 100, -1.5, 1.5);
    fOutputListCharged->Add(fHistPhiEta);

    fhTrigTracks = new TH2D("fhTrigTracks", "fhTrigTracks; pT (trig); PVz", fPtBinsTrigCharged.size() - 1, fPtBinsTrigCharged.data(), 10, -10, 10);
    fOutputListCharged->Add(fhTrigTracks);
    
    //mixing
    fPoolMgr = new AliEventPoolManager(fPoolMaxNEvents, fPoolMinNTracks, fNCentBins,fCentBins.data(), fNzVtxBins, fzVtxBins.data());
    
    if (!fPoolMgr) { AliError("Event Pool manager not created!"); return; }
    fPoolMgr->SetTargetValues(fPoolMinNTracks, 0.1, 5);

    //sparses
    
    Int_t nSteps = 1;
    Double_t binning_deta_tpctpc[33] = {-1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0,    0.1,  0.2,  0.3,  0.4,  0.5, 0.6,  0.7,  0.8,  0.9,  1.0,  1.1,  1.2,  1.3,  1.4,  1.5, 1.6};
    
    Double_t binning_dphi[73] = { -1.570796, -1.483530, -1.396263, -1.308997, -1.221730, -1.134464, -1.047198, -0.959931, -0.872665, -0.785398, -0.698132, -0.610865, -0.523599, -0.436332, -0.349066, -0.261799, -0.174533, -0.087266, 0.0,       0.087266,  0.174533,  0.261799,  0.349066,  0.436332, 0.523599,  0.610865,  0.698132,  0.785398,  0.872665,  0.959931, 1.047198,  1.134464,  1.221730,  1.308997,  1.396263,  1.483530, 1.570796,  1.658063,  1.745329,  1.832596,  1.919862,  2.007129, 2.094395,  2.181662,  2.268928,  2.356194,  2.443461,  2.530727, 2.617994,  2.705260,  2.792527,  2.879793,  2.967060,  3.054326, 3.141593,  3.228859,  3.316126,  3.403392,  3.490659,  3.577925, 3.665191,  3.752458,  3.839724,  3.926991,  4.014257,  4.101524, 4.188790,  4.276057,  4.363323,  4.450590,  4.537856,  4.625123, 4.712389};
    
    const Int_t sizePtTrig = fPtBinsTrigCharged.size() - 1;
    const Int_t sizePtAss = fPtBinsAss.size() - 1;
    const Int_t iBinningTPCTPC[5] = {32,72,sizePtTrig,sizePtAss, 10};
    
    fhChargedSE = new AliTHn("fhChargedSE", "fhChargedSE", nSteps, 5, iBinningTPCTPC);
    
    fhChargedSE->SetBinLimits(0, binning_deta_tpctpc);
    
    fhChargedSE->SetBinLimits(1, binning_dphi);
    
    fhChargedSE->SetBinLimits(2, fPtBinsTrigCharged.data());
    
    fhChargedSE->SetBinLimits(3, fPtBinsAss.data());
    
    fhChargedSE->SetBinLimits(4, -10,10);
    fhChargedSE->SetVarTitle(0, "#Delta#eta");
    fhChargedSE->SetVarTitle(1, "#Delta#phi");
    fhChargedSE->SetVarTitle(2, "p_{T} [GeV/c] (trig)");
    fhChargedSE->SetVarTitle(3, "p_{T} [GeV/c] (ass)");
    fhChargedSE->SetVarTitle(4, "PVz");
    fOutputListCharged->Add(fhChargedSE);
    
    

    fhChargedME = new AliTHn("fhChargedME", "fhChargedME", nSteps, 5, iBinningTPCTPC);
    fhChargedME->SetBinLimits(0, binning_deta_tpctpc);
    fhChargedME->SetBinLimits(1, binning_dphi);
    fhChargedME->SetBinLimits(2, fPtBinsTrigCharged.data());
    fhChargedME->SetBinLimits(3, fPtBinsAss.data());
    fhChargedME->SetBinLimits(4, -10,10);
    fhChargedME->SetVarTitle(0, "#Delta#eta");
    fhChargedME->SetVarTitle(1, "#Delta#phi");
    fhChargedME->SetVarTitle(2, "p_{T} [GeV/c] (trig)");
    fhChargedME->SetVarTitle(3, "p_{T} [GeV/c] (ass)");
    fhChargedME->SetVarTitle(4, "PVz");
    fOutputListCharged->Add(fhChargedME);


    if(fUseEfficiency) {
      fInputListEfficiency = (TList*) GetInputData(1);
      if(fEfficiencyEtaDependent && fAbsEtaMax > 0.8) AliWarning("Efficiency loading -- eta can be out of range!");
    }

    PostData(1, fOutputListCharged);
    
}
//_____________________________________________________________________________
void AliAnalysisFlowTask6::UserExec(Option_t *)
{
    
    fhEventCounter->Fill("Input",1);

    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if(!fAOD) { AliError("Event not loaded."); return; }
    if(!IsEventSelected()) { return; }
    
    fhEventCounter->Fill("Event OK",1);

    Int_t iTracks(fAOD->GetNumberOfTracks());
    if(iTracks < 1 ) {
      AliWarning("No tracks in the event.");
      return;
    }

    //!fTracksTrigCharged = new TObjArray;
    fTracksAss = new TObjArray;
    fTracksTrigCharged = new TObjArray;



    for(Int_t i(0); i < iTracks; i++) {
        AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));
        
        if(!track || !IsTrackSelected(track)) { continue; }
        Double_t trackPt = track->Pt();
        if(trackPt > fPtMinAss && trackPt < fPtMaxAss) fTracksAss->Add((AliAODTrack*)track);
        if(trackPt > fPtMinTrig && trackPt < fPtMaxTrig) {
          fTracksTrigCharged->Add((AliAODTrack*)track);
          fhTrigTracks->Fill(trackPt, fPVz);
            
        }

        //example histogram

        fHistPhiEta->Fill(track->Phi(), track->Eta());
        
    }
    if(!fTracksTrigCharged->IsEmpty()){
        //FillCorrelations();
        //FillCorrelationsMixed();
    }
    
    fTracksTrigCharged->Clear();
        delete fTracksTrigCharged;

    fTracksAss->Clear();
        delete fTracksAss;

    PostData(1, fOutputListCharged);

    
    
}
//_____________________________________________________________________________
void AliAnalysisFlowTask6::Terminate(Option_t *)
{}
//_____________________________________________________________________________
Bool_t AliAnalysisFlowTask6::IsEventSelected()
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

/*
Bool_t AliAnalysisFlowTask6::IsEventSelected()
{
  fhEventCounter->Fill("EventOK",1);

  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) mgr->GetInputEventHandler();
  UInt_t fSelectMask = inputHandler->IsEventSelected();
  if(!(fSelectMask & fTrigger)) { return kFALSE; }
  fhEventCounter->Fill("TriggerOK",1);
    
  if(fIsHMpp) fEventCuts.OverrideAutomaticTriggerSelection(AliVEvent::kHighMultV0, true);
  if(!fEventCuts.AcceptEvent(fAOD)) { return kFALSE; }
  fhEventCounter->Fill("CutsOK",1);

  AliMultSelection* multSelection = (AliMultSelection*) fAOD->FindListObject("MultSelection");
  if(!multSelection) { return kFALSE; }
  fhEventCounter->Fill("MultOK",1);
  Float_t dPercentile = multSelection->GetMultiplicityPercentile(fCentEstimator);
  if(dPercentile > 100 || dPercentile < 0) { return kFALSE; }
  fhEventCounter->Fill("PercOK",1);

  if(fCentMax > 0.0 && (dPercentile < fCentMin || dPercentile > fCentMax)) { return kFALSE; }
  fhEventCounter->Fill("CentOK",1);
  fCentrality = (Double_t) dPercentile;
        

  fPVz = fAOD->GetPrimaryVertex()->GetZ();
  if(TMath::Abs(fPVz) >= 10.0) { return kFALSE; }
  fhEventCounter->Fill("PVzOK",1);

  fbSign = (InputEvent()->GetMagneticField() > 0) ? 1 : -1;

  return kTRUE;
}
 */

//_____________________________________________________________________________
Bool_t AliAnalysisFlowTask6::IsTrackSelected(const AliAODTrack* track) const
{
  if(!track->TestFilterBit(fFilterBit)) { return kFALSE; }
  if(track->GetTPCNcls() < 70 && fFilterBit != 2) { return kFALSE; }
  if(fAbsEtaMax > 0.0 && TMath::Abs(track->Eta()) > fAbsEtaMax) { return kFALSE; }
  if(track->Charge() == 0) { return kFALSE; }

  return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliAnalysisFlowTask6::IsEventRejectedAddPileUp() const
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
/*
Bool_t AliAnalysisFlowTask6::AreWeightsLoaded()
{
  if(!fInputWeights) {AliError("Weights input list not loaded"); return kFALSE; }
  fhWeight = (TH1D*)fInputWeights->FindObject("hWeight");
  if(!fhWeight) {AliError("Weights not loaded"); return kFALSE; }
  return kTRUE;
}
 */
//_____________________________________________________________________________
