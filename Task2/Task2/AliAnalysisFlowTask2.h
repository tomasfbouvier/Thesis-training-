#ifndef AliAnalysisFlowTask2_H
#define AliAnalysisFlowTask2_H

#include "TChain.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TProfile2D.h"
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

class AliAnalysisFlowTask2 : public AliAnalysisTaskSE
{
    public:
                                AliAnalysisFlowTask2();
                                AliAnalysisFlowTask2(const char *name, Bool_t bUseWeights);

        virtual                 ~AliAnalysisFlowTask2();

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
        AliAnalysisFlowTask2(const AliAnalysisFlowTask2&); // not implemented
        AliAnalysisFlowTask2& operator=(const AliAnalysisFlowTask2&); // not implemented

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
        Double_t                fAbsEtaMax; // [1.0]
        Double_t                fPVtxCutZ; // [10.0] PV z cut (in cm)
        Double_t                fPVz; // [99.9] PV z-coordinate
    
        Double_t                fDCAtoVtxZ; // [3.0] max distance to vertex Z
        Double_t                fDCAtoVtxXY; // [3.0] max distance to vertex XY
        UShort_t                fNTPCClusters; // [70] Number of TPC clusters
        UShort_t                fNITSClusters; // [2] Number of ITS clusters
        Double_t                fTPCchi2perClusterMin; // [0.1] Minimum Chi2 per Cluster
        Double_t                fTPCchi2perClusterMax; // [4.] Maximum Chi2 per Cluster

        Double_t                v2;
        Double_t                v3;
        Double_t                v4;
    
    
        TString                 fCentEstimator; // [kV0M]
        AliEventCuts            fEventCuts;

        //outputs
        TList*                  fOutputList;    //! output list
        //output events QA
        TH1D*                   fhEventCounter; //!z
        TH2F*                   fHistPhiEta;    //!
        TProfile*               fMyProfile_v2;     //! is this correct?
        TProfile*               fMyProfile_v3;     //! is this correct?
        TProfile*               fMyProfile_v4;     //! is this correct?
    
        TProfile2D*             fMyProfile_AB;     //! is this correct?
        TProfile2D*             fMyProfile_BA;     //! is this correct?

    
    
        TList*                  fInputWeights;    //! list w weights
        //output events QA
        TH1D*                   fhWeight; //!
    
        Double_t                fgap;

        vector<Double_t>        Q2;
    
        Double_t Total_weight_A;
        Double_t Total_weight_B;
    
    
        struct data{
            vector<vector<Double_t>>    pA;
            vector<Double_t>            QA;
            vector<vector<Double_t>>    pB;
            vector<Double_t>            QB;
            Double_t                    Total_weight_A;
            Double_t                    Total_weight_B;
        } Cumulants;
    
        
        Double_t                weight;
        
        Bool_t                  fWeighting;
        Bool_t                  fGapping;
    
        Int_t                   fNBins;


        void update_Q(vector<Double_t> &Q2, Double_t &Total_weight_A, Double_t &Total_weight_B, Int_t n, Double_t Phi, Double_t Eta, Double_t weight );

        Double_t C2(const vector<Double_t> Q, const vector<Double_t> p);
    
        void update_Cumulants(struct data &Cumulants, const AliAODTrack* track,  const Int_t n, const  Double_t weight);
        
        //!Double_t C2(const vector<Double_t> Q, const vector<Double_t> p);

        void Calculate_corelations(const struct data &Cumulants);

        Int_t CheckPtbin(Double_t pT);
        
    
        ClassDef(AliAnalysisFlowTask2, 1);
};

#endif




//! TODO ask zuzana to provide the file with her results so we can plot them together.
