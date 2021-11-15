#ifdef __CLING__
#include "AliAnalysisAlien.h"
#include "AliAnalysisManager.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisTaskFlowExample2.cxx"
#endif
#include "AliMultSelectionTask.h"


void runFlowExample2()
{
    Bool_t local = 1;
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
    //gInterpreter->LoadMacro("AliAnalysisTaskFlowTest.cxx++g");
    AliAnalysisTaskFlowExample2 *task = reinterpret_cast<AliAnalysisTaskFlowExample2*>(gInterpreter->ExecuteMacro("AddFlowExampleTask2.C(\"FlowTask\", \"weight.root\")"));
#else
    gROOT->LoadMacro("AliAnalysisTaskFlowExample2.cxx++g");
    gROOT->LoadMacro("AddFlowExampleTask2.C");
    AliAnalysisTaskFlowExample2 *task = AddFlowExampleTask2();
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
        chain->Add("AliAOD1.root");
        // chain->Add("AliAOD_pp_LHC16g.root");
        mgr->StartAnalysis("local", chain);
    } else {
        AliAnalysisAlien *alienHandler = new AliAnalysisAlien();
        alienHandler->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");
        alienHandler->SetAdditionalLibs("AliAnalysisTaskFlowExample2.cxx AliAnalysisTaskFlowExample2.h");
        // alienHandler->SetAdditionalLibs("PWGCFFLOWGF.par libPWGEMCALbase.so");
        alienHandler->SetAnalysisSource("AliAnalysisTaskFlowExample2.cxx");
        alienHandler->SetAliPhysicsVersion("vAN-20201214_ROOT6-1");
        // alienHandler->SetAPIVersion("V1.1x");
        alienHandler->SetGridDataDir("/alice/data/2016/LHC16g");
        alienHandler->SetDataPattern("pass1/AOD208/*/AliAOD.root");
        alienHandler->SetRunPrefix("000");
        alienHandler->AddRunNumber(254332);
        alienHandler->SetSplitMaxInputFileNumber(40);
        alienHandler->SetExecutable("FlowExampleTask.sh");
        alienHandler->SetTTL(10000);
        alienHandler->SetJDLName("FlowExampleTask.jdl");

        alienHandler->SetOutputToRunNo(kTRUE);
        alienHandler->SetKeepLogs(kTRUE);
        alienHandler->SetMaxMergeStages(1);
        alienHandler->SetMergeViaJDL(kTRUE);

        alienHandler->SetGridWorkingDir("TEST");
        alienHandler->SetGridOutputDir("myOutputDir");

        mgr->SetGridHandler(alienHandler);
        if(gridTest) {
            alienHandler->SetNtestFiles(1);
            alienHandler->SetRunMode("test");
            mgr->StartAnalysis("grid");
        } else {
            alienHandler->SetRunMode("full");
            mgr->StartAnalysis("grid");
        }
    }
}

