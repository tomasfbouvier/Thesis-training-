#ifndef AliAnalysisFlowTask6_H
#define AliAnalysisFlowTask6_H

#include "TChain.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TList.h"
#include "TMath.h"
#include "TString.h"
#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliMultSelection.h"
#include "AliAODInputHandler.h"
#include "TComplex.h"
#include "AliEventPoolManager.h"
#include "AliTHn.h"

class AliAnalysisFlowTask6 : public AliAnalysisTaskSE
{
    public:
                                AliAnalysisFlowTask6();
                                AliAnalysisFlowTask6(const char *name, Bool_t bUseEff);

        virtual                 ~AliAnalysisFlowTask6();

        virtual void            UserCreateOutputObjects();
        virtual void            UserExec(Option_t* option);
        virtual void            Terminate(Option_t* option);

        //event selection
        void                    SetTrigger(AliVEvent::EOfflineTriggerTypes trigger) { fTrigger = trigger; }
        void                    SetRejectAddPileUp(Bool_t use = kTRUE) { fEventRejectAddPileUp = use; }
        void                    SetCentralityEst(TString est){ fCentEstimator = est; }
        void                    SetCentralityBin(Int_t nOfBins, Double_t cenMin, Double_t cenMax){ fCentBinNum = nOfBins; fCentMin = cenMin; fCentMax = cenMax; }
        //track selection
        void                    SetFilterBit(UInt_t filter) { fFilterBit = filter; }
        void                    SetPtRange(Double_t min, Double_t max) {fPtMin = min; fPtMax = max; }
        void                    SetAbsEta(Double_t etaAbs) {fAbsEtaMax = etaAbs; }

    
    
    private:
        AliAnalysisFlowTask6(const AliAnalysisFlowTask6&); // not implemented
        AliAnalysisFlowTask6& operator=(const AliAnalysisFlowTask6&); // not implemented

        Bool_t                  IsEventSelected();
        Bool_t                  IsEventRejectedAddPileUp() const;
        Bool_t                  IsTrackSelected(const AliAODTrack* track) const;
    
        
        Bool_t                  AreWeightsLoaded();
        Double_t                GetWeight(const Double_t dPhi);

        AliAODEvent*            fAOD;           //! input event
        //event and track selection
        AliVEvent::EOfflineTriggerTypes    fTrigger; // [AliVEvent::kINT7]
        Bool_t                  fIsPP; // [kFALSE]
        Bool_t                  fEventRejectAddPileUp; // [kTRUE]
        Bool_t                  fUseWeights; // [kTrue]
        UInt_t                  fFilterBit; // [96]
        Int_t                   fCentBinNum; // [10]
        Double_t                fCentMin; // [0.0]
        Double_t                fCentMax; // [100.0]
        Double_t                fCentrality; // [1]
        Double_t                fPtMin; // [0.2]
        Double_t                fPtMax; // [5.0]
        Double_t                fAbsEtaMin; // [0.4]
        Double_t                fAbsEtaMax; // [0.8]
        Double_t                fPVtxCutZ; // [10.0] PV z cut (in cm)
        Double_t                fPVz; // [99.9] PV z-coordinate
    
        Double_t                fDCAtoVtxZ; // [3.0] max distance to vertex Z
        Double_t                fDCAtoVtxXY; // [3.0] max distance to vertex XY
        UShort_t                fNTPCClusters; // [70] Number of TPC clusters
        UShort_t                fNITSClusters; // [2] Number of ITS clusters
        Double_t                fTPCchi2perClusterMin; // [0.1] Minimum Chi2 per Cluster
        Double_t                fTPCchi2perClusterMax; // [4.] Maximum Chi2 per Cluster
        

        //!Double_t                value;
    
        TString                 fCentEstimator; // [kV0M]
        AliEventCuts            fEventCuts;

        //output events QA
        TH1D*                   fhEventCounter; //!z  //!
    




    
    //! ADDED FOR CORRELATION TASK
    
        std::vector<Double_t>   fzVtxBins;
        Int_t                   fNzVtxBins; // number of PV z bins
        std::vector<Double_t>   fCentBins;
        TList*                  fOutputListCharged;    //! output list
        TH2D*                   fHistPhiEta; //!
        TH2D*                   fhTrigTracks; //!
        std::vector<Double_t>   fPtBinsTrigCharged;
        AliEventPoolManager*    fPoolMgr;  //!  event pool manager for Event Mixing
        Int_t                   fPoolMinNTracks;   // minimum number of tracks in the pool
        std::vector<Double_t>   fPtBinsAss;
        
    
        AliTHn*                 fhChargedSE; //!
        Int_t                   fPoolMaxNEvents;   // maximum number of events in the pool
        Int_t                   fNCentBins; // number of centrality bins
        AliTHn*                 fhChargedME; //!
        Bool_t                  fUseEfficiency; // [kFALSE]
        TList*                  fInputListEfficiency;    //! input list
        Bool_t                  fEfficiencyEtaDependent; // [kFALSE]
    
    //!
    

  
        ClassDef(AliAnalysisFlowTask6, 1);
};

#endif




//! TODO ask zuzana to provide the file with her results so we can plot them together.
