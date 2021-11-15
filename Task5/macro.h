void calculate_SC(TProfile *prof_results_SC, TProfile *prof_results_NSC, TProfile *prof_4pc_m_n, TProfile *prof_2pc_m,  TProfile *prof_2pc_n, TProfile *prof_2pc_m_gap,  TProfile *prof_2pc_n_gap);
TMultiGraph* represent(TList* fOutputList, TList* fErrorList);
TProfile* load_results(const char * fname, const char * particles);
void macro();



TProfile* Mean_profile_SC_5_2;
TProfile* Mean_profile_SC_5_3;
TProfile* Mean_profile_SC_4_3;
TProfile* Mean_profile_SC_4_2;
TProfile* Mean_profile_SC_3_2;

TProfile* Error_profile_SC_5_2;
TProfile* Error_profile_SC_5_3;
TProfile* Error_profile_SC_4_3;
TProfile* Error_profile_SC_4_2;
TProfile* Error_profile_SC_3_2;

TProfile* Mean_profile_NSC_5_2;
TProfile* Mean_profile_NSC_5_3;
TProfile* Mean_profile_NSC_4_3;
TProfile* Mean_profile_NSC_4_2;
TProfile* Mean_profile_NSC_3_2;

TProfile* Error_profile_NSC_5_2;
TProfile* Error_profile_NSC_5_3;
TProfile* Error_profile_NSC_4_3;
TProfile* Error_profile_NSC_4_2;
TProfile* Error_profile_NSC_3_2;


const char *name_SC_5_2[]= {"4pc_5_2_1", "4pc_5_2_2", "4pc_5_2_3", "4pc_5_2_4","4pc_5_2_5", "4pc_5_2_6", "4pc_5_2_7", "4pc_5_2_8", "4pc_5_2_9" ,"4pc_5_2_10"}; // enough to hold all numbers up to 64-bits TODO: FIX THIS FUCKIN SHIT
const char *name_SC_5_3[]= {"4pc_5_3_1", "4pc_5_3_2", "4pc_5_3_3", "4pc_5_3_4","4pc_5_3_5", "4pc_5_3_6", "4pc_5_3_7", "4pc_5_3_8", "4pc_5_3_9" ,"4pc_5_3_10"}; // enough to hold all numbers up to 64-bits TODO: FIX THIS FUCKIN SHIT
const char *name_SC_4_3[]= {"4pc_4_3_1", "4pc_4_3_2", "4pc_4_3_3", "4pc_4_3_4","4pc_4_3_5", "4pc_4_3_6", "4pc_4_3_7", "4pc_4_3_8", "4pc_4_3_9" ,"4pc_4_3_10"}; // enough to hold all numbers up to 64-bits TODO: FIX  FUCKIN SHIT
const char *name_SC_4_2[]= {"4pc_4_2_1", "4pc_4_2_2", "4pc_4_2_3", "4pc_4_2_4","4pc_4_2_5", "4pc_4_2_6", "4pc_4_2_7", "4pc_4_2_8", "4pc_4_2_9" ,"4pc_4_2_10"}; // enough to hold all numbers up to 64-bits TODO: FIX  FUCKIN SHIT
const char *name_SC_3_2[]= {"4pc_3_2_1", "4pc_3_2_2", "4pc_3_2_3", "4pc_3_2_4","4pc_3_2_5", "4pc_3_2_6", "4pc_3_2_7", "4pc_3_2_8", "4pc_3_2_9" ,"4pc_3_2_10"}; // enough to hold all numbers up to 64-bits TODO: FIX  FUCKIN SHIT

const char *name_2pc_2[]= {"2pc_2_1", "2pc_2_2", "2pc_2_3", "2pc_2_4","2pc_2_5", "2pc_2_6", "2pc_2_7", "2pc_2_8", "2pc_2_9" ,"2pc_2_10"}; // enough to hold all numbers up to 64-bits TODO: FIX  FUCKIN SHIT
const char *name_2pc_3[]= {"2pc_3_1", "2pc_3_2", "2pc_3_3", "2pc_3_4","2pc_3_5", "2pc_3_6", "2pc_3_7", "2pc_3_8", "2pc_3_9" ,"2pc_3_10"}; // enough to hold all numbers up to 64-bits TODO: FIX  FUCKIN SHIT
const char *name_2pc_4[]= {"2pc_4_1", "2pc_4_2", "2pc_4_3", "2pc_4_4","2pc_4_5", "2pc_4_6", "2pc_4_7", "2pc_4_8", "2pc_4_9" ,"2pc_4_10"}; // enough to hold all numbers up to 64-bits TODO: FIX  FUCKIN SHIT
const char *name_2pc_5[]= {"2pc_5_1", "2pc_5_2", "2pc_5_3", "2pc_5_4","2pc_5_5", "2pc_5_6", "2pc_5_7", "2pc_5_8", "2pc_5_9" ,"2pc_5_10"}; // enough to hold all numbers up to 64-bits TODO: FIX  FUCKIN SHIT


const char *name_vn[]= {"2_1", "2_2", "2_3", "2_4","2_5", "2_6", "2_7", "2_8", "2_9" ,"2_10"}; // enough to hold all numbers up to 64-bits TODO: FIX THIS FUCKIN SHIT

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


const char *names[]={"SC(5,2)", "SC(5,3)", "SC(4,3)", "SC(4,2)","SC(3,2)"};

TList* fOutputListSC;
TList* fErrorListSC;

TList* fOutputListNSC;
TList* fErrorListNSC;



