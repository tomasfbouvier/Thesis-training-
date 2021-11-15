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

/* AliAnalysisFlowTask3
 *
 * example flow class:
 * event selection
 * track selection
 * loading particle weights if needed (possible to change run (in)dependent)
 * filling Q vectors (only needed combinations of harmonics & powers! check if this is what you need!)
 * calculate m-particle correlations (using generic algorithm)
 */

#include "AliAnalysisFlowTask3.h"

class AliAnalysisFlowTask3;

using namespace std;

ClassImp(AliAnalysisFlowTask3)

AliAnalysisFlowTask3::AliAnalysisFlowTask3() : AliAnalysisTaskSE(),
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
AliAnalysisFlowTask3::AliAnalysisFlowTask3(const char* name, Bool_t bUseWeights) : AliAnalysisTaskSE(name),
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
    maxHarm(4),
    maxWeightPower(4),
    maxbins(99),
    fFlowVecQpos{},
    fFlowVecPpos{},
    fFlowMatPpos{},
    fFlowVecSpos{}
{
    DefineInput(0, TChain::Class());
    if(bUseWeights) { DefineInput(1, TList::Class()); }
    DefineOutput(1, TList::Class());
}
//_____________________________________________________________________________
AliAnalysisFlowTask3::~AliAnalysisFlowTask3()
{
    if(fOutputList) delete fOutputList;
}
//_____________________________________________________________________________
void AliAnalysisFlowTask3::UserCreateOutputObjects()
{

    // OpenFile(1);

    // creating output lists
    fOutputList = new TList();
    fOutputList->SetOwner(kTRUE);
    
    Double_t bin_edges[] = {0.,5.,10.,20.,30.,40.,50.,60.,70.,80.,90.};
    
    TString sEventCounterLabel[] = {"Input", "Trigger OK", "Event cuts OK", "Event OK"};
    const Int_t iEventCounterBins = sizeof(sEventCounterLabel)/sizeof(sEventCounterLabel[0]);
    fhEventCounter = new TH1D("fhEventCounter","Event Counter",iEventCounterBins,0,iEventCounterBins);
    fOutputList->Add(fhEventCounter);

    fMyProfile_v_2_2= new TProfile2D("Big_2","Profile eta vs value; centrality;  pT ",10, bin_edges, maxbins, 0.2, 5.);


    fOutputList->Add(fMyProfile_v_2_2);
    
    fMyProfile_v_2_4= new TProfile2D("Big_4","Profile eta vs value; centrality;  pT ",10, bin_edges, maxbins, 0.2, 5.);
    
    fOutputList->Add(fMyProfile_v_2_4);

    
    const char *name_2[10]= {"2_1", "2_2", "2_3", "2_4","2_5", "2_6", "2_7", "2_8", "2_9" ,"2_10"}; // enough to hold all numbers up to 64-bits TODO: FIX THIS FUCKIN SHIT
    const char *name_4[10]= {"4_1", "4_2", "4_3", "4_4","4_5", "4_6", "4_7", "4_8", "4_9" ,"4_10"}; // enough to hold all numbers up to 64-bits TODO: FIX THIS FUCKIN SHIT

    for(Int_t i(0); i<10; i++){
        fProfiles_collection2[i]= new TProfile2D(name_2[i],"Profile eta vs value; centrality;  pT ",10, bin_edges, maxbins, 0.2, 5.);
        fOutputList->Add(fProfiles_collection2[i]);
    }
    
    for(Int_t i(0); i<10; i++){
        fProfiles_collection4[i]= new TProfile2D(name_4[i],"Profile eta vs value; centrality;  pT ",10, bin_edges, maxbins, 0.2, 5.);
        fOutputList->Add(fProfiles_collection4[i]);
    }
    
   
    
    
    
    /*/ TRIAL OF SWITHC TO CONTAINER
    
    TObjArray *OAforPt=new TObjArray();
    //No gap:
    OAforPt->Add(new TNamed("MidV22","MidV22"));
    for(Int_t i=0;i<maxbins;i++)
      OAforPt->Add(new TNamed(Form("MidV22_pt_%i",i+1),"MidV22_pTDiff"));
    OAforPt->Add(new TNamed("MidV24","MidV24"));
    for(Int_t i=0;i<maxbins;i++)
      OAforPt->Add(new TNamed(Form("MidV24_pt_%i",i+1),"MidV24_pTDiff"));
    
    
    
    fFC = new AliGFWFlowContainer();
    
    fFC->Initialize(OAforPt,7,13,10)
    
    */// END OF TRIAL
    
    if(fUseWeights) fInputWeights = (TList*) GetInputData(1);
    
    PostData(1, fOutputList);
}
//_____________________________________________________________________________
void AliAnalysisFlowTask3::UserExec(Option_t *)
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

    ResetFlowMatrix(fFlowMatPpos, maxHarm, maxWeightPower, maxbins);
    
    for(Int_t i(0); i < iTracks; i++) {
        
        AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));
        if(!track || !IsTrackSelected(track)) { continue; }
        
        if(fWeighting) weight =  fhWeight->GetBinContent(fhWeight->GetBin(track->Phi()));
        else weight=1.;
        
        FillVectors(track, fFlowVecQpos, fFlowMatPpos);
        
    }

    fiHarm[0]=2;fiHarm[1]=-2;fiHarm[2]=2;fiHarm[3]=-2;fiHarm[4]=2;
    
    for(Int_t bin(0); bin<maxbins; bin++){
        Get_subset(fFlowVecPpos, fFlowMatPpos, maxHarm, maxWeightPower, bin);
        Get_subset(fFlowVecSpos, fFlowMatPpos, maxHarm, maxWeightPower, bin);
        
        Calculate_correlation_2(bin);
        Calculate_correlation_4(bin);
    }
    PostData(1, fOutputList);
    
    
    
}
//_____________________________________________________________________________
void AliAnalysisFlowTask3::Terminate(Option_t *)
{}
//_____________________________________________________________________________
Bool_t AliAnalysisFlowTask3::IsEventSelected()
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
Bool_t AliAnalysisFlowTask3::IsTrackSelected(const AliAODTrack* track) const
{
    if(!track->TestFilterBit(fFilterBit)) { return kFALSE; }
    if(track->GetTPCNcls() < 70 && fFilterBit != 2) { return kFALSE; }
    if(fPtMin > 0 && track->Pt() < fPtMin) { return kFALSE; }
    if(fPtMax > 0 && track->Pt() > fPtMax) { return kFALSE; }
    if(fAbsEtaMax > 0 && TMath::Abs(track->Eta()) > fAbsEtaMax) { return kFALSE; }
    if(track->DCA() > fDCAtoVtxXY) { return kFALSE; }
    if(track->ZAtDCA() > fDCAtoVtxZ) { return kFALSE; }
    if(track->GetTPCNcls() < fNTPCClusters)  { return kFALSE; }
    if(track->GetITSNcls()< fNITSClusters)  { return kFALSE; }
    if(track->GetTPCchi2perCluster() < fTPCchi2perClusterMin || track->GetTPCchi2perCluster() > fTPCchi2perClusterMax)  { return kFALSE; }
    
    
    return kTRUE;
}
//_____________________________________________________________________________
Bool_t AliAnalysisFlowTask3::IsEventRejectedAddPileUp() const
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

Bool_t AliAnalysisFlowTask3::AreWeightsLoaded()
{
  if(!fInputWeights) {AliError("Weights input list not loaded"); return kFALSE; }
  fhWeight = (TH1D*)fInputWeights->FindObject("hWeight");
  if(!fhWeight) {AliError("Weights not loaded"); return kFALSE; }
  return kTRUE;
}
//_____________________________________________________________________________
Double_t AliAnalysisFlowTask3::GetWeight(const Double_t dPhi)
{
  if(!fUseWeights) return 1.0;
  if(!fhWeight) {AliFatal("Weights not loaded"); return -1.0; }
  return fhWeight->GetBinContent(fhWeight->FindFixBin(dPhi));
}

// ============================================================================

TComplex AliAnalysisFlowTask3::Q(const Int_t n, const Int_t p) const
{
 if (n < 0) return TComplex::Conjugate(fFlowVecQpos[-n][p]);
 else return fFlowVecQpos[n][p];
}
// ============================================================================
TComplex AliAnalysisFlowTask3::QGapPos(const Int_t n, const Int_t p) const
{
 if (n < 0) return TComplex::Conjugate(fFlowVecQpos[-n][p]);
 else return fFlowVecQpos[n][p];
}
// ============================================================================
TComplex AliAnalysisFlowTask3::QGapNeg(const Int_t n, const Int_t p) const
{
 if(n < 0) return TComplex::Conjugate(fFlowVecQneg[-n][p]);
 else return fFlowVecQneg[n][p];
}


TComplex AliAnalysisFlowTask3::P(const Int_t n, const Int_t p) const
{
  if(n < 0) return TComplex::Conjugate(fFlowVecPpos[-n][p]);
  else return fFlowVecPpos[n][p];
}
// ============================================================================
TComplex AliAnalysisFlowTask3::PGapPos(const Int_t n, const Int_t p) const
{
  if(n < 0) return TComplex::Conjugate(fFlowVecPpos[-n][p]);
  else return fFlowVecPpos[n][p];
}
// ============================================================================
TComplex AliAnalysisFlowTask3::PGapNeg(const Int_t n, const Int_t p) const
{
  if(n < 0) return TComplex::Conjugate(fFlowVecPneg[-n][p]);
  else return fFlowVecPneg[n][p];
}
// ============================================================================
TComplex AliAnalysisFlowTask3::S(const Int_t n, const Int_t p) const
{
  if(n < 0) return TComplex::Conjugate(fFlowVecSpos[-n][p]);
  else return fFlowVecSpos[n][p];
}

void AliAnalysisFlowTask3::ResetFlowVector(TComplex (&array)[fFlowNumHarmonicsMax][fFlowNumWeightPowersMax], Int_t maxHarm, Int_t maxWeightPower)
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

void AliAnalysisFlowTask3::ResetFlowMatrix(TComplex (&array)[fNBins][fFlowNumHarmonicsMax][fFlowNumWeightPowersMax], Int_t maxHarm, Int_t maxWeightPower, Int_t maxbins)
{
  // Reset RFPs (Q) array values to TComplex(0,0,kFALSE) for given array
  // *************************************************************
  for(Int_t iBin(0); iBin <= maxbins; ++iBin) {
      for(Int_t iHarm(0); iHarm <= maxHarm; ++iHarm) {
        for(Int_t iPower(0); iPower <= maxWeightPower; ++iPower) {
            array[iBin][iHarm][iPower](0.0,0.0);
        }
      }
  }
  return;
}

void AliAnalysisFlowTask3::Calculate_correlation_2(Int_t bin){
    

   //! cout<<cNom<<"\t";
    // important to distinguish between different "setup"

    Double_t value = 0.0;
    Double_t valueNeg = 0.0;

    TComplex cDenom = TwoDiff(0,0);
    TComplex cNom = TwoDiff(fiHarm[0],fiHarm[1]);
    
    Double_t bin_center= ((Double_t)(bin))/maxbins*(fPtMax-fPtMin)+fPtMin;
    
    Bool_t bFillPos = kFALSE;


    if(cDenom.Re() > 0.0) { bFillPos = kTRUE; value = cNom.Re() / cDenom.Re(); }
    if(bFillPos && TMath::Abs(value) > 1.0) { bFillPos = kFALSE; }

    if(value > -10 && fCentrality< 90){
        Int_t random_number = rand() % 10 ;
        if(bFillPos) fMyProfile_v_2_2->Fill(fCentrality, bin_center, value, cDenom.Re());
        if(bFillPos) fProfiles_collection2[random_number]->Fill(fCentrality, bin_center, value, cDenom.Re());
    }
    return ;
    
}


void AliAnalysisFlowTask3::Calculate_correlation_4(Int_t bin){
    

   //! cout<<cNom<<"\t";
    // important to distinguish between different "setup"
    
    Double_t value = 0.0;
    Double_t valueNeg = 0.0;

    TComplex cDenom = FourDiff(0,0,0,0);
    TComplex cNom = FourDiff(fiHarm[0],fiHarm[1], fiHarm[2], fiHarm[3]);
    
    Double_t bin_center= ((Double_t)(bin))/maxbins*(fPtMax-fPtMin)+fPtMin;
    
    Bool_t bFillPos = kFALSE;


    if(cDenom.Re() > 0.0) { bFillPos = kTRUE; value = cNom.Re() / cDenom.Re(); }
    if(bFillPos && TMath::Abs(value) > 1.0) { bFillPos = kFALSE; }

    if(value > -10 && fCentrality< 90){
        Int_t random_number = rand() % 10 ;
        if(bFillPos) fMyProfile_v_2_4->Fill(fCentrality, bin_center, value, cDenom.Re());
        if(bFillPos) fProfiles_collection4[random_number]->Fill(fCentrality, bin_center, value, cDenom.Re());
    }
    return ;
    
}

TComplex AliAnalysisFlowTask3::Two(const Int_t n1, const Int_t n2) const
{
 TComplex formula = Q(n1,1)*Q(n2,1) - Q(n1+n2,2);
 return formula;
}

// ============================================================================
/*TComplex AliAnalysisFlowTask3::TwoGap(const Int_t n1, const Int_t n2) const
{
 TComplex formula = QGapPos(n1,1)*QGapNeg(n2,1);
 return formula;
}*/
// ============================================================================
TComplex AliAnalysisFlowTask3::TwoDiff(const Int_t n1, const Int_t n2) const
{
 TComplex formula = P(n1,1)*Q(n2,1) - S(n1+n2,2);
 return formula;
}
// ============================================================================
/*TComplex AliAnalysisFlowTask3::TwoDiffGapPos(const Int_t n1, const Int_t n2) const
{
 TComplex formula = PGapPos(n1,1)*QGapNeg(n2,1);
 return formula;
}
// ============================================================================
TComplex AliAnalysisFlowTask3::TwoDiffGapNeg(const Int_t n1, const Int_t n2) const
{
 TComplex formula = PGapNeg(n1,1)*QGapPos(n2,1);
 return formula;
}
*/
TComplex AliAnalysisFlowTask3::Four(const Int_t n1, const Int_t n2, const Int_t n3, const Int_t n4) const
{
  TComplex formula = Q(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)-Q(n1+n2,2)*Q(n3,1)*Q(n4,1)-Q(n2,1)*Q(n1+n3,2)*Q(n4,1)
                    - Q(n1,1)*Q(n2+n3,2)*Q(n4,1)+2.0*Q(n1+n2+n3,3)*Q(n4,1)-Q(n2,1)*Q(n3,1)*Q(n1+n4,2)
                    + Q(n2+n3,2)*Q(n1+n4,2)-Q(n1,1)*Q(n3,1)*Q(n2+n4,2)+Q(n1+n3,2)*Q(n2+n4,2)
                    + 2.0*Q(n3,1)*Q(n1+n2+n4,3)-Q(n1,1)*Q(n2,1)*Q(n3+n4,2)+Q(n1+n2,2)*Q(n3+n4,2)
                    + 2.0*Q(n2,1)*Q(n1+n3+n4,3)+2.0*Q(n1,1)*Q(n2+n3+n4,3)-6.0*Q(n1+n2+n3+n4,4);
  return formula;
}



void AliAnalysisFlowTask3::FillVectors(const AliAODTrack* track, TComplex (&fFlowVecQpos)[fFlowNumHarmonicsMax][fFlowNumWeightPowersMax], TComplex (&fFlowMatPpos)[fNBins][fFlowNumHarmonicsMax][fFlowNumWeightPowersMax]){
    
    //!                         I asume that from now on I will always work with gap

    Int_t Ptbin = CheckPtbin(track->Pt());
      
        // RFP in positive eta acceptance

    for(Int_t iHarm(0); iHarm <= maxHarm; iHarm++){
        for(Int_t iPower(0); iPower <= maxWeightPower; iPower++){
            Double_t dCos = TMath::Power(weight,iPower) * TMath::Cos(iHarm * track->Phi());
            Double_t dSin = TMath::Power(weight,iPower) * TMath::Sin(iHarm * track->Phi());
            fFlowVecQpos[iHarm][iPower] += TComplex(dCos,dSin,kFALSE);
            fFlowMatPpos[Ptbin][iHarm][iPower] += TComplex(dCos,dSin,kFALSE);
        }
    }

    return ;
    
};




Int_t AliAnalysisFlowTask3::CheckPtbin(Double_t Pt)
{
    Int_t bin= (Int_t) ((Pt-fPtMin)/(fPtMax-fPtMin)*maxbins);
   //! cout<<Pt<<"\t"<<bin<<"\n";
    return bin;
};
void AliAnalysisFlowTask3::Get_subset (TComplex (&array)[fFlowNumHarmonicsMax][fFlowNumWeightPowersMax], TComplex (&array2)[fNBins][fFlowNumHarmonicsMax][fFlowNumWeightPowersMax], Int_t maxHarm, Int_t maxWeightPower, Int_t bin)
{
  // Reset RFPs (Q) array values to TComplex(0,0,kFALSE) for given array
  // *************************************************************
  for(Int_t iHarm(0); iHarm <= maxHarm; ++iHarm) {
    for(Int_t iPower(0); iPower <= maxWeightPower; ++iPower) {
      array[iHarm][iPower]= array2[bin][iHarm][iPower];
    }
  }
  return;
}

TComplex AliAnalysisFlowTask3::FourDiff(const Int_t n1, const Int_t n2, const Int_t n3, const Int_t n4) const
{
  TComplex formula = P(n1,1)*Q(n2,1)*Q(n3,1)*Q(n4,1)-S(n1+n2,2)*Q(n3,1)*Q(n4,1)-Q(n2,1)*S(n1+n3,2)*Q(n4,1)
                    - P(n1,1)*Q(n2+n3,2)*Q(n4,1)+2.0*S(n1+n2+n3,3)*Q(n4,1)-Q(n2,1)*Q(n3,1)*S(n1+n4,2)
                    + Q(n2+n3,2)*S(n1+n4,2)-P(n1,1)*Q(n3,1)*Q(n2+n4,2)+S(n1+n3,2)*Q(n2+n4,2)
                    + 2.0*Q(n3,1)*S(n1+n2+n4,3)-P(n1,1)*Q(n2,1)*Q(n3+n4,2)+S(n1+n2,2)*Q(n3+n4,2)
                    + 2.0*Q(n2,1)*S(n1+n3+n4,3)+2.0*P(n1,1)*Q(n2+n3+n4,3)-6.0*S(n1+n2+n3+n4,4);
  return formula;
}
