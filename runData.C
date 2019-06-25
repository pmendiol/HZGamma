{
	gSystem->SetBuildDir("tmpdir", kTRUE);
	gROOT->LoadMacro("/afs/cern.ch/work/p/pmendiol/2018Yields/HtoZG-master/MYcode/plugin/RoccoR.cc++");
	gROOT->LoadMacro("Analysis.C+");
	system("mkdir -p minitree");

	int RegType = 2; 
	
	const char* inpathsMu[] = {
	"/data1/ggNtuples/V08_00_26_05/job_DoubleMu_Run2016B_FebReminiAOD/ggtree_data_*.root",
	"/data1/ggNtuples/V08_00_26_05/job_DoubleMu_Run2016C_FebReminiAOD/ggtree_data_*.root",
	"/data1/ggNtuples/V08_00_26_05/job_DoubleMu_Run2016D_FebReminiAOD/ggtree_data_*.root",
	"/data1/ggNtuples/V08_00_26_05/job_DoubleMu_Run2016E_FebReminiAOD/ggtree_data_*.root",
	"/data1/ggNtuples/V08_00_26_05/job_DoubleMu_Run2016F_FebReminiAOD*/ggtree_data_*.root",
	"/data1/ggNtuples/V08_00_26_05/job_DoubleMu_Run2016G_FebReminiAOD/ggtree_data_*.root",
	"/data1/ggNtuples/V08_00_26_05/job_DoubleMu_Run2016H_FebReminiAODv*/ggtree_data_*.root",
	};
	
	const char* inpathsEle[] = {
	"/data1/ggNtuples/V08_00_26_05/job_DoubleEG_Run2016B_FebReminiAOD/ggtree_data_*.root",
	"/data1/ggNtuples/V08_00_26_05/job_DoubleEG_Run2016C_FebReminiAOD/ggtree_data_*.root",
	"/data1/ggNtuples/V08_00_26_05/job_DoubleEG_Run2016D_FebReminiAOD/ggtree_data_*.root",
	"/data1/ggNtuples/V08_00_26_05/job_DoubleEG_Run2016E_FebReminiAOD/ggtree_data_*.root",
	"/data1/ggNtuples/V08_00_26_05/job_DoubleEG_Run2016F_FebReminiAOD*/ggtree_data_*.root",
	"/data1/ggNtuples/V08_00_26_05/job_DoubleEG_Run2016G_FebReminiAOD/ggtree_data_*.root",
	"/data1/ggNtuples/V08_00_26_05/job_DoubleEG_Run2016H_FebReminiAODv*/ggtree_data_*.root",
	};

	char outname[500];
	const char* data_period[7]={"B","C","D","E","F","G","H"};


	for(int n=0;n<7;n++){
		sprintf(outname,"minitree/MuMu2016_Data%s.root",data_period[n]);
		xAna(inpathsMu[n], outname, 1,0.205814347,35.9,0,RegType,true);
	} 
	
	cout << "MuMu with Regression is finished" << endl;
	cout << "==========================================================" << endl;
	for(int n=0;n<7;n++){
		sprintf(outname,"minitree/EE2016_Data%s.root",data_period[n]);
		xAna(inpathsEle[n], outname, 0,0.205814347,35.9,0,RegType,true);
	}
	cout << "EE with Regression is finished" << endl;


	
	
	
	
}

