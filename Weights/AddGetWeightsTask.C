AliAnalysisGetWeights* AddGetWeightsTask(TString name = "FlowTask", const char* suffix = "")
{
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        return 0x0;
    }
    if (!mgr->GetInputEventHandler()) {
        return 0x0;
    }

    TString fileName = AliAnalysisManager::GetCommonFileName();
    fileName += ":FlowExampleTask";
    TString taskName = Form("%s_%s",name.Data(),suffix);

    AliAnalysisGetWeights* task = new AliAnalysisGetWeights(taskName.Data());
    if(!task) return 0x0;

    task->SelectCollisionCandidates(AliVEvent::kAnyINT);
    mgr->AddTask(task);
    mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
    mgr->ConnectOutput(task,1,mgr->CreateContainer("Output", TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));

    return task;
}
