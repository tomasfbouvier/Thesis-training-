#ifdef __CLING__
#include "AliAnalysisAlien.h"
#include "AliAnalysisManager.h"
#include "AliAODInputHandler.h"
#include "AliMultSelectionTask.h"
#include "AliAnalysisFlowTask2.h"
#endif
#include "AliMultSelectionTask.h"


void runFlowTask2()
{
    Bool_t local = 0;
    Bool_t gridTest = 0;

#if !defined (__CINT__) || defined (__CLING__)
    gInterpreter->ProcessLine(".include $ROOTSYS/include");
    gInterpreter->ProcessLine(".include $ALICE_ROOT/include");
#else
    gROOT->ProcessLine(".include $ROOTSYS/include");
    gROOT->ProcessLine(".include $ALICE_ROOT/include");
#endif

    AliAnalysisManager *mgr = new AliAnalysisManager("AnalysisTaskExample");
    AliAODInputHandler *aodH = new AliAODInputHandler();
    mgr->SetInputEventHandler(aodH);

    AliMultSelectionTask* taskMultSelection = reinterpret_cast<AliMultSelectionTask*>(gInterpreter->ProcessLine(Form(".x %s(kFALSE)", gSystem->ExpandPathName("$ALICE_PHYSICS/OADB/COMMON/MULTIPLICITY/macros/AddTaskMultSelection.C"))));
    taskMultSelection->SetSelectedTriggerClass(AliVEvent::kINT7);
    // taskMultSelection->SetSelectedTriggerClass(AliVEvent::kHighMultV0);

#if !defined (__CINT__) || defined (__CLING__)
    gInterpreter->LoadMacro("AliAnalysisFlowTask2.cxx++g");
    //!AliAnalysisFlowTask2 *task = reinterpret_cast<AliAnalysisFlowTask2*>(gInterpreter->ExecuteMacro("AddFlowTask2.C(\"FlowTask\", \"/Users/tomasfernandezbouvier/Desktop/Thesis-training-/Weights/weight.root\")"));
    AliAnalysisFlowTask2 *task = reinterpret_cast<AliAnalysisFlowTask2*>(gInterpreter->ExecuteMacro("AddFlowTask2.C(\"FlowTask\", \"alien:///alice/cern.ch/user/f/fernandt/weight.root\")"));

#else
    gROOT->LoadMacro("AliAnalysisFlowTask2.cxx++g");
    gROOT->LoadMacro("AddFlowTask2.C");
    AliAnalysisFlowTask2 *task = AddFlowTask2();
#endif

    // task->SetTrigger(AliVEvent::kHighMultV0);
    task->SetTrigger(AliVEvent::kINT7);
    // task->AddFlowTask({2,-2});
    // task->AddFlowTask({2,-2});
    // task->AddFlowTask({2,2,-2,-2});
    // task->SetMaxPowersVector({5,0,4,0,3});
    // task->SetNofSamples(5);

    if(!mgr->InitAnalysis()) return;
    mgr->SetDebugLevel(2);
    mgr->PrintStatus();
    mgr->SetUseProgressBar(1, 25);

    if(local) {
        TChain* chain = new TChain("aodTree");
        chain->Add("AliAOD.root");
        // chain->Add("AliAOD_pp_LHC16g.root");
        mgr->StartAnalysis("local", chain);
    }
    else {
        AliAnalysisAlien *alienHandler = new AliAnalysisAlien();
        alienHandler->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");
        alienHandler->SetAdditionalLibs("AliAnalysisFlowTask2.cxx AliAnalysisFlowTask2.h");
        // alienHandler->SetAdditionalLibs("PWGCFFLOWGF.par libPWGEMCALbase.so");
        alienHandler->SetAnalysisSource("AliAnalysisFlowTask2.cxx");
        // select the aliphysics version. all other packages
        // are LOADED AUTOMATICALLY!
        alienHandler->SetAliPhysicsVersion("vAN-20210906_ROOT6-1");
        // set the Alien API version
        alienHandler->SetAPIVersion("V1.1x");
        // select the input data
        alienHandler->SetGridDataDir("/alice/data/2015/LHC15o");
        alienHandler->SetDataPattern("*pass1/AOD194/*/AliAOD.root");
        // MC has no prefix, data has prefix 000
        alienHandler->SetRunPrefix("000");
        // runnumber
        alienHandler->AddRunNumber(244918);
        // number of files per subjob
        alienHandler->SetSplitMaxInputFileNumber(40);
        alienHandler->SetExecutable("myTask.sh");
        // specify how many seconds your job may take
        alienHandler->SetTTL(50000);
        alienHandler->SetJDLName("myTask.jdl");

        alienHandler->SetOutputToRunNo(kTRUE);
        alienHandler->SetKeepLogs(kTRUE);
        // merging: run with kTRUE to merge on grid
        // after re-running the jobs in SetRunMode("terminate")
        // (see below) mode, set SetMergeViaJDL(kFALSE)
        // to collect final results
        alienHandler->SetMaxMergeStages(1);
        alienHandler->SetMergeViaJDL(kTRUE);

        // define the output folders
        alienHandler->SetGridWorkingDir("Task2");
        alienHandler->SetGridOutputDir("Task2OutputDir");

        // connect the alien plugin to the manager
        mgr->SetGridHandler(alienHandler);
        if(gridTest) {
            // speficy on how many files you want to run
            alienHandler->SetNtestFiles(1);
            // and launch the analysis
            alienHandler->SetRunMode("test");
            mgr->StartAnalysis("grid");
        } else {
            // else launch the full grid analysis
            alienHandler->SetRunMode("fuñll");
            mgr->StartAnalysis("grid");
        }
    }
}

