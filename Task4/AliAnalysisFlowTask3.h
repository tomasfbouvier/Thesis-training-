#ifndef AliAnalysisFlowTask3_H
#define AliAnalysisFlowTask3_H

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
#include "TComplex.h"
#include "AliGFWFlowContainer.h"

class AliAnalysisFlowTask3 : public AliAnalysisTaskSE
{
    public:
                                AliAnalysisFlowTask3();
                                AliAnalysisFlowTask3(const char *name, Bool_t bUseWeights);

        virtual                 ~AliAnalysisFlowTask3();

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
        AliAnalysisFlowTask3(const AliAnalysisFlowTask3&); // not implemented
        AliAnalysisFlowTask3& operator=(const AliAnalysisFlowTask3&); // not implemented

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
        TProfile2D*               fMyProfile_v_2_2;     //! i
    
        TProfile2D*               fMyProfile_v_2;   //! is this correct?
        TProfile2D*               fMyProfile_v_2_4;     //! is this correct?
    
    
        //!
    //!
    //!
        TProfile2D*             fProfiles_collection2[10];     //! is this correct?
        TProfile2D*             fProfiles_collection4[10];     //! is this correct?
        
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
        Int_t                   maxbins;
        
        static const Int_t      fFlowNumWeightPowersMax = 5; // maximum weight power length of flow vector array
        static const Int_t      fFlowNumHarmonicsMax = 5;
        static const Int_t      fNBins = 100;
        Int_t                   fiHarm[fFlowNumHarmonicsMax];
    
        TComplex                fFlowVecQpos[fFlowNumHarmonicsMax][fFlowNumWeightPowersMax]; // flow vector array for flow calculation
        TComplex                fFlowVecQneg[fFlowNumHarmonicsMax][fFlowNumWeightPowersMax]; // flow vector array for flow calculation
        TComplex                fFlowVecPpos[fFlowNumHarmonicsMax][fFlowNumWeightPowersMax]; // flow vector array for flow calculation
        TComplex                fFlowVecPneg[fFlowNumHarmonicsMax][fFlowNumWeightPowersMax]; // flow vector array for flow calculation
    
        TComplex                fFlowVecSpos[fFlowNumHarmonicsMax][fFlowNumWeightPowersMax]; // flow vector array for flow calculation
        
        TComplex                fFlowVecSneg[fFlowNumHarmonicsMax][fFlowNumWeightPowersMax]; // flow vector array for flow calculation
        TComplex                fFlowMatPpos[fNBins][fFlowNumHarmonicsMax][fFlowNumWeightPowersMax]; // flow matrix array for flow calculation
        
            
        
        TComplex Two(const Int_t n1, const Int_t n2) const;
        TComplex TwoGap(const Int_t n1, const Int_t n2) const;
        TComplex TwoDiff(const Int_t n1, const Int_t n2) const;
        TComplex TwoDiffGapPos(const Int_t n1, const Int_t n2) const;
        TComplex TwoDiffGapNeg(const Int_t n1, const Int_t n2) const;
        TComplex Q(const Int_t n, const Int_t p) const;
        TComplex QGapPos(const Int_t n, const Int_t p) const;
        TComplex QGapNeg(const Int_t n, const Int_t p) const;
        
        TComplex P(const Int_t n, const Int_t p) const;
        TComplex PGapPos(const Int_t n, const Int_t p) const;
        TComplex PGapNeg(const Int_t n, const Int_t p) const;
    
        TComplex S(const Int_t n, const Int_t p) const;
    

    
        TComplex Four(const Int_t n1, const Int_t n2, const Int_t n3, const Int_t n4) const;
        TComplex Six(const Int_t n1, const Int_t n2, const Int_t n3, const Int_t n4, const Int_t n5, const Int_t n6) const;
        TComplex FourDiff(const Int_t n1, const Int_t n2, const Int_t n3, const Int_t n4) const;
        
        void Calculate_correlation_2(Int_t bin);
        void Calculate_correlation_4(Int_t bin);
            
        
        Int_t   CheckPtbin(Double_t Pt);
        void    FillVectors(const AliAODTrack* track, TComplex (&fFlowVecQpos)[fFlowNumHarmonicsMax][fFlowNumWeightPowersMax], TComplex (&fFlowMatPpos)[fNBins][fFlowNumHarmonicsMax][fFlowNumWeightPowersMax]);

        void    ResetFlowVector(TComplex (&array)[fFlowNumHarmonicsMax][fFlowNumWeightPowersMax], Int_t maxHarm, Int_t maxWeightPower); // set values to TComplex(0,0,0) for given array

        void    ResetFlowMatrix(TComplex (&array)[fNBins][fFlowNumHarmonicsMax][fFlowNumWeightPowersMax], Int_t maxHarm, Int_t maxWeightPower, Int_t maxbins);
        void    Get_subset (TComplex (&array)[fFlowNumHarmonicsMax][fFlowNumWeightPowersMax], TComplex (&array2)[fNBins][fFlowNumHarmonicsMax][fFlowNumWeightPowersMax], Int_t maxHarm, Int_t maxWeightPower, Int_t bin);

        AliGFWFlowContainer     fFC;
    
        ClassDef(AliAnalysisFlowTask3, 1);
};

#endif
