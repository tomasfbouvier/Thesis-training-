void calculate_SC(TProfile *prof_results_SC, TProfile *prof_results_NSC, TProfile *prof_4pc_m_n, TProfile *prof_2pc_m,  TProfile *prof_2pc_n, TProfile *prof_2pc_m_gap,  TProfile *prof_2pc_n_gap);
void represent(TList* fOutputList, TList* fErrorList);
TProfile* load_results(const char * fname, const char * particles);
void macro();




const char *name_vn[]= {"2_1", "2_2", "2_3", "2_4","2_5", "2_6", "2_7", "2_8", "2_9" ,"2_10", "results_value_2"}; // enough to hold all numbers up to 64-bits TODO: FIX THIS FUCKIN SHIT

/*
TFile* inputFile;
TFile* inputFile2;
TFile* inputFile3;
TFile* inputFile4;
TFile* inputFile5;
TFile* inputFile6;

TList* inputList
TList* inputList2
TList* inputList3
TList* inputList4
TList* inputList5
TList* inputList6
*/


const char *names[]={"SC(5,2)", "SC(5,3)", "SC(4,3)", "SC(4,2)","SC(3,2)", "asdfjj", "asdkooie"};

TList* fMeanList_rho=new TList();
TList* fErrorList_rho=new TList();

TList* fMeanList_chi=new TList();
TList* fErrorList_chi=new TList();

static std::string profile_nom_names[] = {"prof_3pc_4_22","prof_3pc_5_23","prof_4pc_6_222", "prof_3pc_6_24","prof_3pc_6_33", "prof_4pc_7_233", "prof_4pc_8_233"};
static std::string profile_denom_names[] = {"prof_4pc_2", "prof_4pc_23", "prof_6pc_2", "prof_4pc_24", "prof_4pc_3", "prof_6pc_3_2222", "prof_6pc_2_3333"};

const char *profile_nom_names2[] = {"chi4_22;V0M (%); chi","chi5_23;V0M (%); chi","chi6_222;V0M (%); chi", "chi6_24;V0M (%); chi","chi6_33;V0M (%); chi", "chi7_233;V0M (%); chi", "chi8_233;V0M (%); chi"};
static const char *inputFile_names[]={"~/Desktop/Thesis-training-/Task1/results_v4/AnalysisResultsGapWeights.root", "~/Desktop/Thesis-training-/Task1/results_v5/AnalysisResultsGapWeights.root","~/Desktop/Thesis-training-/Task1/results_v6/AnalysisResultsGapWeights.root","~/Desktop/Thesis-training-/Task1/results_v6/AnalysisResultsGapWeights.root","~/Desktop/Thesis-training-/Task1/results_v6/AnalysisResultsGapWeights.root","~/Desktop/Thesis-training-/Task1/results_v7/AnalysisResultsGapWeights.root", "~/Desktop/Thesis-training-/Task1/results_v8/AnalysisResultsGapWeights.root"};
