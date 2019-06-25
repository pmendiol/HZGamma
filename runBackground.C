{
	gSystem->SetBuildDir("tmpdir", kTRUE);
	gROOT->LoadMacro("/afs/cern.ch/work/p/pmendiol/2018Yields/HtoZG-master/MYcode/plugin/RoccoR.cc++");
	gROOT->LoadMacro("Analysis.C+");
	system("mkdir -p minitree");

	int RegType = 2;  
	

	// /*BACKGROUND*/

	/*TTbar*/
	/*MUON CHANNEL*/ 
	xAna("/data6/V08_00_26_09/job_summer16_TT_powheg/*/*/*/000*/ggtree_mc_*.root", 
	"minitree/TTbar_muon.root",1,831760,35.9,0,RegType,true);
	cout << "MuMu with Regression is finished TTbar" << endl;
	cout << "-----------------------------------------" << endl;
	
	/*ELECTRON CHANNEL*/ 
	xAna("/data6/V08_00_26_09/job_summer16_TT_powheg/*/*/*/000*/ggtree_mc_*.root",
	"minitree/TTbar_electron.root",0,831760,35.9,0,RegType,true);
	cout << "MuMu with Regression is finished TTbar" << endl;
	cout << "-----------------------------------------" << endl;
	
	/*SM ZG aMC@NLO*/
	/*MUON CHANNEL*/ 
	xAna("/data6/ggNtuples/V08_00_26_05/job_summer16_Zg_aMCatNLO/ggtree_mc_*.root", 
	"minitree/ZgaMCatNLO_muon.root",1,117864,35.9,1,RegType,true);
	cout << "MuMu with Regression is finished Zg_aMCatNLO" << endl;
	cout << "-----------------------------------------" << endl;
	
	/*ELECTRON CHANNEL*/ 
	xAna("/data6/ggNtuples/V08_00_26_05/job_summer16_Zg_aMCatNLO/ggtree_mc_*.root",
	"minitree/ZgaMCatNLO_electron.root",0,117864,35.9,1,RegType,true);
	cout << "MuMu with Regression is finished Zg_aMCatNLO" << endl;
	cout << "-----------------------------------------" << endl;
	
	/*WZTo2L2Q*/
	/*MUON CHANNEL*/ 
	xAna("/data6/ggNtuples/V08_00_26_05/job_summer16_WZTo2L2Q/ggtree_mc_*.root", 
	"minitree/WZTo2L2Q_muon.root",1,5595,35.9,0,RegType,true);
	cout << "MuMu with Regression is finished WZTo2L2Q" << endl;
	cout << "-----------------------------------------" << endl;
	
	/*ELECTRON CHANNEL*/ 
	xAna("/data6/ggNtuples/V08_00_26_05/job_summer16_WZTo2L2Q/ggtree_mc_*.root",
	"minitree/WZTo2L2Q_electron.root",0,5595,35.9,0,RegType,true);
	cout << "MuMu with Regression is finished WZTo2L2Q" << endl;
	cout << "-----------------------------------------" << endl;
	
	/*WWToLNuQQ*/
	/*MUON CHANNEL*/ 
	xAna("/data6/ggNtuples/V08_00_26_05/job_summer16_WWToLNuQQ/ggtree_mc_*.root", 
	"minitree/WWToLNuQQ_muon.root",1,49997,35.9,0,RegType,true);
	cout << "MuMu with Regression is finished WWToLNuQQ" << endl;
	cout << "-----------------------------------------" << endl;
	
	/*ELECTRON CHANNEL*/ 
	xAna("/data6/ggNtuples/V08_00_26_05/job_summer16_WWToLNuQQ/ggtree_mc_*.root",
	"minitree/WWToLNuQQ_electron.root",0,49997,35.9,0,RegType,true);
	cout << "MuMu with Regression is finished WWToLNuQQ" << endl;
	cout << "-----------------------------------------" << endl;
	
	/*ZZTo2L2Q*/
	/*MUON CHANNEL*/ 
	xAna("/data6/ggNtuples/V08_00_26_05/job_summer16_ZZTo2L2Q/ggtree_mc_*.root", 
	"minitree/ZZTo2L2Q_muon.root",1,3220,35.9,0,RegType,true);
	cout << "MuMu with Regression is finished ZZTo2L2Q" << endl;
	cout << "-----------------------------------------" << endl;
	
	/*ELECTRON CHANNEL*/ 
	xAna("/data6/ggNtuples/V08_00_26_05/job_summer16_ZZTo2L2Q/ggtree_mc_*.root",
	"minitree/ZZTo2L2Q_electron.root",0,3220,35.9,0,RegType,true);
	cout << "MuMu with Regression is finished ZZTo2L2Q" << endl;
	cout << "-----------------------------------------" << endl;
	
	/*DY+JETS*/
	/*MUON CHANNEL*/ 
	xAna("/data6/ggNtuples/V08_00_26_05/job_summer16_DYJetsToLL_m50_aMCatNLO/ggtree_mc_*.root", 
	"minitree/DYJets_muon.root",1,5943200,35.9,1,RegType,true);
	cout << "MuMu with Regression is finished DY+Jets" << endl;
	cout << "-----------------------------------------" << endl;
	
	/*ELECTRON CHANNEL*/ 
	xAna("/data6/ggNtuples/V08_00_26_05/job_summer16_DYJetsToLL_m50_aMCatNLO/ggtree_mc_*.root", 
	"minitree/DYJets_electron.root",0,5943200,35.9,1,RegType,true);
	cout << "MuMu with Regression is finished DY+Jets" << endl;
	cout << "-----------------------------------------" << endl;
	
	
	
	
}