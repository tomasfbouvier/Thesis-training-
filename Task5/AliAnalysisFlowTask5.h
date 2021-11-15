#ifndef AliAnalysisFlowTask5_H
#define AliAnalysisFlowTask5_H

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

class AliAnalysisFlowTask5 : public AliAnalysisTaskSE
{
    public:
                                AliAnalysisFlowTask5();
                                AliAnalysisFlowTask5(const char *name, Bool_t bUseWeights);

        virtual                 ~AliAnalysisFlowTask5();

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
        AliAnalysisFlowTask5(const AliAnalysisFlowTask5&); // not implemented
        AliAnalysisFlowTask5& operator=(const AliAnalysisFlowTask5&); // not implemented

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
        

        //!Double_t                value;
    
        TString                 fCentEstimator; // [kV0M]
        AliEventCuts            fEventCuts;

        //outputs
        TList*                  fOutputList;    //! output list
        //output events QA
        TH1D*                   fhEventCounter; //!z  //!
    
        TProfile*             fProfiles_collection_SC_5_2[11];     //!
        TProfile*             fProfiles_collection_SC_5_3[11];     //!
        TProfile*             fProfiles_collection_SC_4_3[11];     //!
        TProfile*             fProfiles_collection_SC_4_2[11];     //!
        TProfile*             fProfiles_collection_SC_3_2[11];     //!
    
        TProfile*             fProfiles_collection_2pc_2[11]; //!
        TProfile*             fProfiles_collection_2pc_3[11]; //!
        TProfile*             fProfiles_collection_2pc_4[11]; //!
        TProfile*             fProfiles_collection_2pc_5[11]; //!
        
    


        TList*                  fInputWeights;    //! list w weights
        //output events QA
        TH1D*                   fhWeight; //!
    
        Double_t                fgap;

        vector<Double_t>        Q2;
    
        Double_t Total_weight_A;
        Double_t Total_weight_B;
    
        Double_t                weight;
        
        Bool_t                  fWeighting;
        Bool_t                  fGapping;
        Bool_t                  fDiff;
    
        Int_t                   maxHarm;
        Int_t                   maxWeightPower;
    
    

            
        static const Int_t      fFlowNumWeightPowersMax = 6; // maximum weight power length of flow vector array
        static const Int_t      fFlowNumHarmonicsMax = 6; // maximum harmonics length of flow vector array
    
        Int_t                   fiHarm[fFlowNumHarmonicsMax];
        TComplex                fFlowVecQpos[fFlowNumHarmonicsMax][fFlowNumWeightPowersMax]; // flow vector array for flow calculation
        TComplex                fFlowVecQneg[fFlowNumHarmonicsMax][fFlowNumWeightPowersMax]; // flow vector array for flow calculation

        
        
        TComplex TwoGap(const Int_t n1, const Int_t n2) const;
        TComplex Q(const Int_t n, const Int_t p) const;
        TComplex QGapPos(const Int_t n, const Int_t p) const;
        TComplex QGapNeg(const Int_t n, const Int_t p) const;

        TComplex Two(const Int_t n1, const Int_t n2) const;
        TComplex Four(const Int_t n1, const Int_t n2, const Int_t n3, const Int_t n4) const;

        void Calculate_correlation_2(TProfile*  filling_collect[]);
        void Calculate_correlation_4(TProfile*  filling_collect[]);

        void    ResetFlowVector(TComplex (&array)[fFlowNumHarmonicsMax][fFlowNumWeightPowersMax], Int_t maxHarm, Int_t maxWeightPower); // set values to TComplex(0,0,0) for given array

        
    
        ClassDef(AliAnalysisFlowTask5, 1);
};

#endif




//! TODO ask zuzana to provide the file with her results so we can plot them together.
