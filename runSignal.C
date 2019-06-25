{
	
	gSystem->SetBuildDir("tmpdir", kTRUE);
	gROOT->LoadMacro("/afs/cern.ch/work/p/pmendiol/2018Yields/HtoZG-master/MYcode/plugin/RoccoR.cc++");
	gROOT->LoadMacro("Analysis.C+");
	system("mkdir -p minitree");

	int RegType = 2;


	/*MUON CHANNEL*/ 
	xAna("/data6/ggNtuples/V08_00_26_05/job_summer16_HZg_ZH_120GeV/ggtree_mc_*.root",
	"minitree/MuMu2016_Regressed_MC_ZH_120GeV.root",1,7.79E-03,35.9,0,RegType,true);
	cout << "MuMu with Regression is finished ZH 120GeV" << endl;
	cout << "-----------------------------------------" << endl;
	
	xAna("/data6/ggNtuples/V08_00_26_05/job_summer16_HZg_ZH_125GeV/ggtree_mc_*.root", 
	"minitree/MuMu2016_Regressed_MC_ZH_125GeV.root",1,6.93E-03,35.9,0,RegType,true);
	cout << "MuMu with Regression is finished ZH 125GeV" << endl;
	cout << "-----------------------------------------" << endl;
	
	xAna("/data6/ggNtuples/V08_00_26_05/job_summer16_HZg_ZH_130GeV/ggtree_mc_*.root", 
	"minitree/MuMu2016_Regressed_MC_ZH_130GeV.root",1,6.19E-03,35.9,0,RegType,true);
	cout << "MuMu with Regression is finished ZH 130GeV" << endl;
	cout << "-----------------------------------------" << endl;


	cout << "===================================================" << endl;
	
	// /*ELECTRON CHANNEL*/
	xAna("/data6/ggNtuples/V08_00_26_05/job_summer16_HZg_ZH_120GeV/ggtree_mc_*.root", 
	"minitree/EE2016_Regressed_MC_ZH_120GeV.root",0,7.78E-03,35.9,0,RegType,true);
	cout << "EE with Regression is finished ZH 120GeV" << endl;
	cout << "-----------------------------------------" << endl;
	
	xAna("/data6/ggNtuples/V08_00_26_05/job_summer16_HZg_ZH_125GeV/ggtree_mc_*.root", 
	"minitree/EE2016_Regressed_MC_ZH_125GeV.root",0,6.92E-03,35.9,0,RegType,true);
	cout << "EE with Regression is finished ZH 125GeV" << endl;
	cout << "-----------------------------------------" << endl;

	xAna("/data6/ggNtuples/V08_00_26_05/job_summer16_HZg_ZH_130GeV/ggtree_mc_*.root", 
	"minitree/EE2016_Regressed_MC_ZH_130GeV.root",0,6.19E-03,35.9,0,RegType,true);
	cout << "EE with Regression is finished ZH 130GeV" << endl;
	cout << "-----------------------------------------" << endl;
}

// XS*BranchingRatios

// MUON CHANNEL

// 120
// ZH: 7.79E-03
// W+H: 2.37E-02
// W-H: 1.51E-02
// ttH: 2.38E-03

// 125
// ZH: 6.93E-03
// W+H: 2.08E-02
// W-H: 1.32E-02
// ttH: 2.12E-03

// 130
// ZH: 6.19E-03
// W+H: 1.84E-02
// W-H: 1.16E-02
// ttH: 1.90E-03


// ELECTRON CHANNEL

// 120
// ZH: 7.78E-03
// W+H: 2.38E-02
// W-H: 1.52E-02
// ttH: 2.35E-03

// 125
// ZH: 6.92E-03
// W+H: 2.09E-02
// W-H: 1.33E-02
// ttH: 2.09E-03

// 130
// ZH: 6.19E-03
// W+H: 1.85E-02
// W-H: 1.17E-02
// ttH: 1.87E-03

//MUON CHANNEL
// xAna("/data6/ggNtuples/V08_00_26_05/job_summer16_HZg_WplusH_130GeV/ggtree_mc_*.root", 
// "minitree/MuMu2016_Regressed_MC_WplusH.root",1,2.08E-02,35.9,0,RegType,true);
// cout << "MuMu with Regression is finished WplusH" << endl;
// cout << "-----------------------------------------" << endl;

// xAna("/data6/ggNtuples/V08_00_26_05/job_summer16_HZg_WminusH_130GeV/ggtree_mc_*.root", 
// "minitree/MuMu2016_Regressed_MC_WminusH.root",1,1.32E-02,35.9,0,RegType,true);
// cout << "MuMu with Regression is finished WminusH" << endl;
// cout << "-----------------------------------------" << endl;

// xAna("/data6/ggNtuples/V08_00_26_05/job_summer16_HZg_ttH_130GeV/ggtree_mc_*.root", 
// "minitree/MuMu2016_Regressed_MC_ttH.root",1,2.12E-03,35.9,0,RegType,true);
// cout << "MuMu with Regression is finished ttH" << endl; 

//ELECTRON CHANNEL
// xAna("/data6/ggNtuples/V08_00_26_05/job_summer16_HZg_WplusH_130GeV/ggtree_mc_*.root",
// "minitree/EE2016_Regressed_MC_WplusH.root",0,2.09E-02,35.9,0,RegType,true);
// cout << "EE with Regression is finished WplusH" << endl;
// cout << "-----------------------------------------" << endl;

// xAna("/data6/ggNtuples/V08_00_26_05/job_summer16_HZg_WminusH_130GeV/ggtree_mc_*.root",
// "minitree/EE2016_Regressed_MC_WminusH.root",0,1.33E-02,35.9,0,RegType,true);
// cout << "EE with Regression is finished WminusH" << endl;  
// cout << "-----------------------------------------" << endl;

// xAna("/data6/ggNtuples/V08_00_26_05/job_summer16_HZg_ttH_130GeV/ggtree_mc_*.root", 
// "minitree/EE2016_Regressed_MC_ttH.root",0,2.09E-03,35.9,0,RegType,true);
// cout << "EE with Regression is finished ttH" << endl; 
