#include "TFile.h"
#include "FV.h"
#include "TVector2.h"
#include "TVector3.h"
#include <iostream>
#include "TLine.h"
#include "TTree.h"
#include <TCanvas.h>
#include <TMVA/Reader.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>



void multipion_analysis(){            //first bracket

  const int ncuts=10;
  double Selected[ncuts][6][10]={0};

  TString CutsName[ncuts];
  CutsName[0]="EventsInTrueFV";
  CutsName[1]="InFiducialVol";
  CutsName[2]="Topological";
  CutsName[3]="MuonCandidate";
  CutsName[4]="ContainedPion";
  CutsName[5]="MuonIn3Planes";
  CutsName[6]="PionIn3Planes";
  CutsName[7]="ShowerCut";
  CutsName[8]="OpeningAngle";
  CutsName[9]="Nonprotons";


  const int nvariables=26;
  TString Variable[nvariables];
  Variable[0]="TopologicalScore";
  Variable[1]="MuonTrackLength";
  Variable[2]="PionTrackLength";
  Variable[3]="MuPiOpeningangle";
  Variable[4]="MuonPID";
  Variable[5]="PionPID";
  Variable[6]="MuonVtxDistance";
  Variable[7]="PionVtxDistance";
  Variable[8]="MuonPionDistance";
  Variable[9]="MuonUPlaneHits";
  Variable[10]="MuonYPlaneHits";
  Variable[11]="MuonVPlaneHits";
  Variable[12]="PionUPlaneHits";
  Variable[13]="PionYPlaneHits";
  Variable[14]="PionVPlaneHits";
  Variable[15]="ShowerUPlaneHits";
  Variable[16]="ShowerYPlaneHits";
  Variable[17]="ShowerVPlaneHits";
  Variable[18]="NPrimaryShowers";
  Variable[19]="NPrimaryTracks";
  Variable[20]="ShowerVtxDistance";
  Variable[21]="PionTMVAMip";
  Variable[22]="PionTMVAPi";
  Variable[23]="MuonTMVAMip";
  Variable[24]="PionTrkEnergy";
  Variable[25]="MuonTrkEnrgy";

  const int nsamples =6;
  TString Sample[nsamples];
  Sample[0]="Signal";
  Sample[1]="Background";
  Sample[2]="EXT";
  Sample[3]="OOFV";
  Sample[4]="Data";
  Sample[5]="MC";



  const int npions=10;
  TString Pions[npions];
  Pions[0]="0 pions";
  Pions[1]="1 pions";
  Pions[2]="2 pions";
  Pions[3]="3 pions";
  Pions[4]="4 pions";
  Pions[5]="5 pions";
  Pions[6]="6 pions";
  Pions[7]="7 pions";
  Pions[8]="8 pions";
  Pions[9]="9 pions";


  TH1F *Histos[ncuts][nsamples][nvariables][npions];

  for(int c=0;c<ncuts;c++){
    for(int v=0;v<nvariables;v++){
      for(int s=0;s<nsamples;s++){
	 for(int p=0;p<npions; p++){
	
        Histos[c][s][v][p] = new TH1F (Variable[v]+"_"+CutsName[c]+"_"+Sample[s]+"_"+Pions[p],"",100,-1,-1);

	if(s!=4) Histos[c][s][v][p] ->SetLineColor(s+1);
	Histos[c][s][v][p] ->GetXaxis()->SetTitle(Variable[v]);
	Histos[c][s][v][p] ->SetTitle(CutsName[c]);

	 }
      }
    }
  }





  std::vector<std::string> Files = {"/data/uboone/new_numi_flux/Run1_fhc_new_numi_flux_fhc_pandora_ntuple.root" // //MC overlay - new numi flux
   				   ,"/data/uboone/new_numi_flux/Run2_fhc_new_numi_flux_fhc_pandora_ntuple.root"
			           ,"/data/uboone/new_numi_flux/Run4_fhc_new_numi_flux_fhc_pandora_ntuple.root"
  			           ,"/data/uboone/new_numi_flux/Run5_fhc_new_numi_flux_fhc_pandora_ntuple.root"
			           ,"/data/uboone/EXT/neutrinoselection_filt_run1_beamoff.root" // EXT beam off
 		                   ,"/data/uboone/dirt/prodgenie_numi_uboone_overlay_dirt_fhc_mcc9_run1_v28_all_snapshot.root" // Dirt - the old numi flux
			           //,"/data/uboone/beam_on/neutrinoselection_filt_run1_beamon_beamgood.root" // Data -
  };






  std::ofstream output("backtracked_showers.txt");


  //POTs and Trigger from Patrick
  //
  // Run 1 FHC: 2.192e+20 POT, 5748692.0 triggers
  // Run 1 FHC: 3.283e+20 POT, 9846635.0 triggers  // with open trigger
  // Run 2 FHC: 1.268e+20 POT, 3535129.0 triggers
  // Run 4 FHC: 2.075e+20 POT, 4131149.0 triggers
  // Run 5 FHC: 2.231e+20 POT, 5154196.0 triggers

  //Beam Off:
  // Run 1 FHC Triggers 4582248.27
  // Run 2 FHC:         10021551.40
  // Run 4 FHC:         8060024.70 (4c) + 17432345.70 (4d)
  // Run 5 FHC:         19256341.475


  double dataPoT[4] ={0};   // in e+20 PoT
  dataPoT[0] = 2.192+3.283; //Run1 inclusind open trigger
  dataPoT[1] = 1.268;  //Run2
  dataPoT[2] = 2.075;  //Run4
  dataPoT[3] = 2.231;  //Run5

  double MCPoT[4] ={0};     // ine+20
  MCPoT[0] = 23.282;  //Run1
  MCPoT[1] = 24.9337; //Run2
  MCPoT[2] = 11.04+17.295;//Run4
  MCPoT[3] = 63.145; //Run5


  double totalDataPoT=0;
  double totalMCPoT=0;

  double Scale[7] ={1};
  for(int i=0;i<4;i++){
    Scale[i] = dataPoT[i]/MCPoT[i];      // data/simulation;
    totalDataPoT+=dataPoT[i];
    totalMCPoT+=MCPoT[i];
  }


  Scale[5] = totalDataPoT/16.739;      // Run1 dirt from David 1.67392e+21
  Scale[4]  = 5748692.0/9199232.0;  // Run1 beam off triggers scaling to beam on triggers


  double test_min=0.0;
  double test_max=0.03;
  double test_step=0.1;

  for(double test=test_min; test<test_max; test+=test_step)
    { // TEST

      for(int cc=0;cc<ncuts;cc++){
	for(int j=0;j<7;j++){
           for(int p=0;p<npions;p++){
    		  Selected[cc][j][p]=0;
	   }
	}
      }




      for(size_t i_f=0;i_f<Files.size();i_f++){
	TFile* f = TFile::Open(Files.at(i_f).c_str());
	if ((!f) || f->IsZombie()) { delete f; return; } // just a precaution
	TTree *t;

	f->GetObject("nuselection/NeutrinoSelectionFilter", t);


	std::cout << "Processing file: " << Files.at(i_f) << std::endl;
	if (!t) {
	  delete f;
	  std::cout << "No tree NeutrinoSelectionFilter found in file " << Files.at(i_f) << std::endl;
	  return;
	} // just a precaution

	t->SetMakeClass(1); // tree in MakeClass mode
	t->SetBranchStatus("*", 0); // disable all branches

	Int_t           run;
	Int_t           sub;
	Int_t           evt;
	Int_t           nu_pdg;
	Int_t           ccnc;
	Int_t           interaction;
	Float_t         true_nu_px;
	Float_t         true_nu_py;
	Float_t         true_nu_pz;
	Float_t         reco_nu_vtx_sce_x;
	Float_t         reco_nu_vtx_sce_y;
	Float_t         reco_nu_vtx_sce_z;
	vector<int>     *backtracked_pdg;
	vector<float>   *backtracked_start_x;
	vector<float>   *backtracked_start_y;
	vector<float>   *backtracked_start_z;
	Int_t           n_tracks;
	Int_t           n_showers;
	vector<float>   *trk_score_v;
	vector<int>     *mc_pdg;
	vector<float>   *mc_E;
	vector<float>   *mc_vx;
	vector<float>   *mc_vy;
	vector<float>   *mc_vz;
	vector<float>   *mc_endx;
	vector<float>   *mc_endy;
	vector<float>   *mc_endz;
	vector<float>   *mc_px;
	vector<float>   *mc_py;
	vector<float>   *mc_pz;
	vector<float>   *trk_dir_x_v;
	vector<float>   *trk_dir_y_v;
	vector<float>   *trk_dir_z_v;
	vector<float>   *trk_sce_start_x_v;
	vector<float>   *trk_sce_start_y_v;
	vector<float>   *trk_sce_start_z_v;
	vector<float>   *trk_sce_end_x_v;
	vector<float>   *trk_sce_end_y_v;
	vector<float>   *trk_sce_end_z_v;
	vector<float>   *trk_len_v;
	vector<float>   *trk_calo_energy_u_v;
	vector<float>   *trk_calo_energy_v_v;
	vector<float>   *trk_calo_energy_y_v;
	vector<float>   *trk_llr_pid_score_v;
	Float_t         topological_score;
	vector<int>     *pfnplanehits_U;
	vector<int>     *pfnplanehits_V;
	vector<int>     *pfnplanehits_Y;
	Size_t          *trk_id;
	vector<unsigned int> *pfp_generation_v;
	vector<float>   *trk_bragg_p_v;
	vector<float>   *trk_bragg_mu_v;
	vector<float>   *trk_bragg_mip_v;
	vector<float>   *shr_start_x_v;
	vector<float>   *shr_start_y_v;
	vector<float>   *shr_start_z_v;
	vector<float>   *shr_dist_v;
	vector<float>   *shr_moliere_avg_v;
	vector<float>   *shr_moliere_rms_v;
	vector<float>   *shr_tkfit_dedx_u_v;
	vector<float>   *shr_tkfit_dedx_v_v;
	vector<float>   *shr_tkfit_dedx_y_v;





	float trk_bragg_p_v_tmva;
	float trk_bragg_mu_v_tmva;
	float trk_bragg_mip_v_tmva;
	float trk_llr_pid_score_v_tmva;
	float trk_len_v_tmva;
	float trk_sce_end_x_v_tmva;
	float trk_sce_end_y_v_tmva;
	float trk_sce_end_z_v_tmva;
	float trk_score_v_tmva;
	// float isContained_tmva;

	float trk_bragg_p_v_tmva_mu;
	float trk_bragg_mu_v_tmva_mu;
	float trk_bragg_mip_v_tmva_mu;
	float trk_llr_pid_score_v_tmva_mu;
	float trk_len_v_tmva_mu;
	float trk_sce_end_x_v_tmva_mu;
	float trk_sce_end_y_v_tmva_mu;
	float trk_sce_end_z_v_tmva_mu;
	float trk_score_v_tmva_mu;

	float trk_bragg_p_v_tmva_pi;
	float trk_bragg_mu_v_tmva_pi;
	float trk_bragg_mip_v_tmva_pi;
	float trk_llr_pid_score_v_tmva_pi;
	float trk_len_v_tmva_pi;
	float trk_sce_end_x_v_tmva_pi;
	float trk_sce_end_y_v_tmva_pi;
	float trk_sce_end_z_v_tmva_pi;
	float trk_score_v_tmva_pi;


	TMVA::Reader * tmvaReader = new TMVA::Reader();
	tmvaReader->AddVariable("trk_bragg_p_v", &trk_bragg_p_v_tmva);
	tmvaReader->AddVariable("trk_bragg_mu_v", &trk_bragg_mu_v_tmva);
	tmvaReader->AddVariable("trk_bragg_mip_v", &trk_bragg_mip_v_tmva);
	tmvaReader->AddVariable("trk_llr_pid_score_v", &trk_llr_pid_score_v_tmva);
	tmvaReader->AddVariable("trk_score_v", &trk_score_v_tmva);
	tmvaReader->AddVariable("trk_len_v", &trk_len_v_tmva);
	tmvaReader->AddVariable("trk_sce_end_x_v", &trk_sce_end_x_v_tmva);
	tmvaReader->AddVariable("trk_sce_end_y_v", &trk_sce_end_y_v_tmva);
	tmvaReader->AddVariable("trk_sce_end_z_v", &trk_sce_end_z_v_tmva);
	//tmvaReader->AddVariable("isContained_out", &isContained_tmva);

	TMVA::Reader * tmvaReader_mu = new TMVA::Reader();
	tmvaReader_mu->AddVariable("trk_bragg_p_v", &trk_bragg_p_v_tmva_mu);
	tmvaReader_mu->AddVariable("trk_bragg_mu_v", &trk_bragg_mu_v_tmva_mu);
	tmvaReader_mu->AddVariable("trk_bragg_mip_v", &trk_bragg_mip_v_tmva_mu);
	tmvaReader_mu->AddVariable("trk_llr_pid_score_v", &trk_llr_pid_score_v_tmva_mu);
	tmvaReader_mu->AddVariable("trk_score_v", &trk_score_v_tmva_mu);
	tmvaReader_mu->AddVariable("trk_len_v", &trk_len_v_tmva_mu);
	tmvaReader_mu->AddVariable("trk_sce_end_x_v", &trk_sce_end_x_v_tmva_mu);
	tmvaReader_mu->AddVariable("trk_sce_end_y_v", &trk_sce_end_y_v_tmva_mu);
	tmvaReader_mu->AddVariable("trk_sce_end_z_v", &trk_sce_end_z_v_tmva_mu);

	TMVA::Reader * tmvaReader_pi = new TMVA::Reader();
	tmvaReader_pi->AddVariable("trk_bragg_p_v", &trk_bragg_p_v_tmva_pi);
	tmvaReader_pi->AddVariable("trk_bragg_mu_v", &trk_bragg_mu_v_tmva_pi);
	tmvaReader_pi->AddVariable("trk_bragg_mip_v", &trk_bragg_mip_v_tmva_pi);
	tmvaReader_pi->AddVariable("trk_llr_pid_score_v", &trk_llr_pid_score_v_tmva_pi);
	tmvaReader_pi->AddVariable("trk_score_v", &trk_score_v_tmva_pi);
	tmvaReader_pi->AddVariable("trk_len_v", &trk_len_v_tmva_pi);
	tmvaReader_pi->AddVariable("trk_sce_end_x_v", &trk_sce_end_x_v_tmva_pi);
	tmvaReader_pi->AddVariable("trk_sce_end_y_v", &trk_sce_end_y_v_tmva_pi);
	tmvaReader_pi->AddVariable("trk_sce_end_z_v", &trk_sce_end_z_v_tmva_pi);




	std::string weightFileName = "booster_decision_tree/dataset_MIP_BDT/weights/TMVAClassification_BDT.weights.xml";
	tmvaReader->BookMVA("BDT",  weightFileName.c_str());

	std::string weightFileName_mu = "booster_decision_tree/dataset_muon_BDT/weights/TMVAClassification_BDT.weights.xml";
	tmvaReader_mu->BookMVA("BDT",  weightFileName_mu.c_str());

	std::string weightFileName_pi = "booster_decision_tree/dataset_pion_BDT/weights/TMVAClassification_BDT.weights.xml";
	tmvaReader_pi->BookMVA("BDT",  weightFileName_pi.c_str());









	backtracked_pdg = 0;
	backtracked_start_x = 0;
	backtracked_start_y = 0;
	backtracked_start_z = 0;
	trk_score_v = 0;
	mc_pdg = 0;
	mc_E = 0;
	mc_vx = 0;
	mc_vy = 0;
	mc_vz = 0;
	mc_endx = 0;
	mc_endy = 0;
	mc_endz = 0;
	mc_px = 0;
	mc_py = 0;
	mc_pz = 0;
	trk_sce_start_x_v = 0;
	trk_sce_start_y_v = 0;
	trk_sce_start_z_v = 0;
	trk_sce_end_x_v = 0;
	trk_sce_end_y_v = 0;
	trk_sce_end_z_v = 0;
	trk_dir_x_v = 0;
	trk_dir_y_v = 0;
	trk_dir_z_v = 0;
	trk_len_v = 0;
	trk_calo_energy_u_v = 0;
	trk_calo_energy_v_v = 0;
	trk_calo_energy_y_v = 0;
	trk_llr_pid_score_v = 0;
	pfnplanehits_U = 0;
	pfnplanehits_V = 0;
	pfnplanehits_Y = 0;
	trk_id = 0;
	pfp_generation_v = 0;
	trk_bragg_p_v = 0;
	trk_bragg_mu_v = 0;
	trk_bragg_mip_v = 0;
	shr_start_x_v = 0;
	shr_start_y_v = 0;
	shr_start_z_v = 0;
	shr_dist_v = 0;
	shr_moliere_avg_v = 0;
	shr_moliere_rms_v = 0;
	shr_tkfit_dedx_u_v = 0;
	shr_tkfit_dedx_v_v = 0;
	shr_tkfit_dedx_y_v = 0;





	t->SetBranchStatus("run",1);
	t->SetBranchStatus("sub",1);
	t->SetBranchStatus("evt",1);
	t->SetBranchStatus("nu_pdg",1);
	t->SetBranchStatus("ccnc",1);
	t->SetBranchStatus("interaction",1);
	t->SetBranchStatus("true_nu_px",1);
	t->SetBranchStatus("true_nu_py",1);
	t->SetBranchStatus("true_nu_pz",1);
	t->SetBranchStatus("reco_nu_vtx_sce_x",1);
	t->SetBranchStatus("reco_nu_vtx_sce_y",1);
	t->SetBranchStatus("reco_nu_vtx_sce_z",1);
	t->SetBranchStatus("backtracked_pdg");
	t->SetBranchStatus("backtracked_start_x");
	t->SetBranchStatus("backtracked_start_y");
	t->SetBranchStatus("backtracked_start_z");
	t->SetBranchStatus("n_tracks",1);
	t->SetBranchStatus("n_showers",1);
	t->SetBranchStatus("trk_score_v",1);
	t->SetBranchStatus("mc_pdg",1);
	t->SetBranchStatus("mc_E",1);
	t->SetBranchStatus("mc_vx",1);
	t->SetBranchStatus("mc_vy",1);
	t->SetBranchStatus("mc_vz",1);
	t->SetBranchStatus("mc_endx",1);
	t->SetBranchStatus("mc_endy",1);
	t->SetBranchStatus("mc_endz",1);
	t->SetBranchStatus("mc_px",1);
	t->SetBranchStatus("mc_py",1);
	t->SetBranchStatus("mc_pz",1);
	t->SetBranchStatus("trk_sce_start_x_v",1);
	t->SetBranchStatus("trk_sce_start_y_v",1);
	t->SetBranchStatus("trk_sce_start_z_v",1);
	t->SetBranchStatus("trk_sce_end_x_v",1);
	t->SetBranchStatus("trk_sce_end_y_v",1);
	t->SetBranchStatus("trk_sce_end_z_v",1);
	t->SetBranchStatus("trk_dir_x_v",1);
	t->SetBranchStatus("trk_dir_y_v",1);
	t->SetBranchStatus("trk_dir_z_v",1);
	t->SetBranchStatus("trk_len_v",1);
	t->SetBranchStatus("trk_calo_energy_u_v",1);
	t->SetBranchStatus("trk_calo_energy_v_v",1);
	t->SetBranchStatus("trk_calo_energy_y_v",1);
	t->SetBranchStatus("trk_llr_pid_score_v",1);
	t->SetBranchStatus("topological_score",1);
	t->SetBranchStatus("pfnplanehits_U",1);
	t->SetBranchStatus("pfnplanehits_V",1);
	t->SetBranchStatus("pfnplanehits_Y",1);
	t->SetBranchStatus("trk_id",1);
	t->SetBranchStatus("pfp_generation_v",1);
	t->SetBranchStatus("trk_bragg_p_v",1);
	t->SetBranchStatus("trk_bragg_mu_v",1);
	t->SetBranchStatus("trk_bragg_mip_v",1);
	t->SetBranchStatus("shr_start_x_v",1);
	t->SetBranchStatus("shr_start_y_v",1);
	t->SetBranchStatus("shr_start_z_v",1);
	t->SetBranchStatus("shr_dist_v",1);
	t->SetBranchStatus("backtracked_purity",1);
	t->SetBranchStatus("shr_moliere_avg_v",1);
	t->SetBranchStatus("shr_moliere_rms_v",1);
	t->SetBranchStatus("shr_tkfit_dedx_u_v",1);
	t->SetBranchStatus("shr_tkfit_dedx_v_v",1);
	t->SetBranchStatus("shr_tkfit_dedx_y_v",1);







	t->SetBranchAddress("run", &run);
	t->SetBranchAddress("sub", &sub);
	t->SetBranchAddress("evt", &evt);
	t->SetBranchAddress("nu_pdg", &nu_pdg);
	t->SetBranchAddress("ccnc", &ccnc);
	t->SetBranchAddress("interaction", &interaction);
	t->SetBranchAddress("true_nu_px", &true_nu_px);
	t->SetBranchAddress("true_nu_py", &true_nu_py);
	t->SetBranchAddress("true_nu_pz", &true_nu_pz);
	t->SetBranchAddress("reco_nu_vtx_sce_x", &reco_nu_vtx_sce_x);
	t->SetBranchAddress("reco_nu_vtx_sce_y", &reco_nu_vtx_sce_y);
	t->SetBranchAddress("reco_nu_vtx_sce_z", &reco_nu_vtx_sce_z);
	t->SetBranchAddress("backtracked_pdg",&backtracked_pdg);
	t->SetBranchAddress("backtracked_start_x",&backtracked_start_x);
	t->SetBranchAddress("backtracked_start_y",&backtracked_start_y);
	t->SetBranchAddress("backtracked_start_z",&backtracked_start_z);
	t->SetBranchAddress("n_tracks", &n_tracks);
	t->SetBranchAddress("n_showers", &n_showers);
	t->SetBranchAddress("trk_score_v", &trk_score_v);
	t->SetBranchAddress("mc_pdg", &mc_pdg);
	t->SetBranchAddress("mc_E", &mc_E);
	t->SetBranchAddress("mc_vx", &mc_vx);
	t->SetBranchAddress("mc_vy", &mc_vy);
	t->SetBranchAddress("mc_vz", &mc_vz);
	t->SetBranchAddress("mc_endx", &mc_endx);
	t->SetBranchAddress("mc_endy", &mc_endy);
	t->SetBranchAddress("mc_endz", &mc_endz);
	t->SetBranchAddress("mc_px", &mc_px);
	t->SetBranchAddress("mc_py", &mc_py);
	t->SetBranchAddress("mc_pz", &mc_pz);
	t->SetBranchAddress("trk_sce_start_x_v", &trk_sce_start_x_v);
	t->SetBranchAddress("trk_sce_start_y_v", &trk_sce_start_y_v);
	t->SetBranchAddress("trk_sce_start_z_v", &trk_sce_start_z_v);
	t->SetBranchAddress("trk_sce_end_x_v", &trk_sce_end_x_v);
	t->SetBranchAddress("trk_sce_end_y_v", &trk_sce_end_y_v);
	t->SetBranchAddress("trk_sce_end_z_v", &trk_sce_end_z_v);
	t->SetBranchAddress("trk_dir_x_v", &trk_dir_x_v);
	t->SetBranchAddress("trk_dir_y_v", &trk_dir_y_v);
	t->SetBranchAddress("trk_dir_z_v", &trk_dir_z_v);
	t->SetBranchAddress("trk_len_v", &trk_len_v);
	t->SetBranchAddress("trk_calo_energy_u_v", &trk_calo_energy_u_v);
	t->SetBranchAddress("trk_calo_energy_v_v", &trk_calo_energy_v_v);
	t->SetBranchAddress("trk_calo_energy_y_v", &trk_calo_energy_y_v);
	t->SetBranchAddress("trk_llr_pid_score_v", &trk_llr_pid_score_v);
	t->SetBranchAddress("topological_score", &topological_score);
	t->SetBranchAddress("pfnplanehits_U", &pfnplanehits_U);
	t->SetBranchAddress("pfnplanehits_V", &pfnplanehits_V);
	t->SetBranchAddress("pfnplanehits_Y", &pfnplanehits_Y);
	t->SetBranchAddress("trk_id", &trk_id);
	t->SetBranchAddress("pfp_generation_v", &pfp_generation_v);
	t->SetBranchAddress("trk_bragg_p_v", &trk_bragg_p_v);
	t->SetBranchAddress("trk_bragg_mu_v", &trk_bragg_mu_v);
	t->SetBranchAddress("trk_bragg_mip_v", &trk_bragg_mip_v);
	t->SetBranchAddress("shr_start_x_v", &shr_start_x_v);
	t->SetBranchAddress("shr_start_y_v", &shr_start_y_v);
	t->SetBranchAddress("shr_start_z_v", &shr_start_z_v);
	t->SetBranchAddress("shr_dist_v", &shr_dist_v);
	t->SetBranchAddress("shr_moliere_avg_v", &shr_moliere_avg_v);
	t->SetBranchAddress("shr_moliere_rms_v", &shr_moliere_rms_v);
	t->SetBranchAddress("shr_tkfit_dedx_u_v",&shr_tkfit_dedx_u_v);
	t->SetBranchAddress("shr_tkfit_dedx_v_v",&shr_tkfit_dedx_v_v);
	t->SetBranchAddress("shr_tkfit_dedx_y_v",&shr_tkfit_dedx_y_v);




	const Long64_t nentries =t->GetEntries();
	int tracknumber =0;


	std::cout<<"Number of Events="<<nentries<<std::endl;
	for(int ientry=0;ientry<nentries;ientry++){
	  t->GetEntry(ientry);

	  if(ientry%100000==0)std::cout<<ientry<<std::endl;

	  int pion_number = 0;
	  int muon_index = -1;
	  int pion_index = -1;
	  int shower_index = -1;
	  bool sel_has_muon_candidate = false;
	  bool sel_has_pion_candidate = false;
	  bool sel_nu_mu_cc = false;
	  bool sel_passed_topo_cut = false;
	  double muon_trk_score = 0.8;
	  double muon_trk_len = 10.0; //cm
	  double muon_pid_score = 0.2;
	  double topo_cut = 0.9;//0.67;
	  bool pion_in_gap = false;
	  bool muon_in_gap = false;
	  int contained_pions = 0;
	  bool sel_contained_pions = false;
	  double mu_pi_opening_angle = 0.0;
	  //double muon_is_correct = 0.0;
	  //double sel_has_a_muon = 0.0;
	  bool one_shower = false;
	  bool plane_hits = false;
	  bool mol_avg = false;
	  bool shower_cut = false;
	  bool Cuts[10][10] ={false};



	  //      if(mc_pdg->size() == 0) continue;



	  bool signal[10] = {false};

	  bool in_fiducial_volume_true = false;

	  int nprotons=0, nchargedpions=0, nmuons=0, npionszero=0, npionsmin=0;
	  if(mc_pdg->size() ) {


	    for(size_t i_mc=0; i_mc<mc_pdg->size();i_mc++){

	      if(mc_pdg->at(i_mc) == 2212)     nprotons++;
	      if(abs(mc_pdg->at(i_mc)) == 13)  nmuons++;
	      if(abs(mc_pdg->at(i_mc)) == 211) nchargedpions++;
	      if(mc_pdg->at(i_mc) == 111)      npionszero++;
	    }


	    TVector3 Reco_Vertex(reco_nu_vtx_sce_x,reco_nu_vtx_sce_y,reco_nu_vtx_sce_z);
	    TVector3 true_vertex(mc_vx->at(0),mc_vy->at(0),mc_vz->at(0));

	    in_fiducial_volume_true = inFV(true_vertex);
	    bool in_fiducial_volume_reco = inFV(Reco_Vertex);

	    //Signal definition
	    //      bool signal = false;


	    if( (in_fiducial_volume_true) && (nmuons==1) && (npionszero==0)  && (abs(nu_pdg) == 14) ){
//              if(nchargedpions>0) std::cout<<nchargedpions<<std::endl;
	      signal[nchargedpions] = true; //
	    }

	  }




	  //reconstructed staff
	  //
	  //



	  TVector3 Reco_Vertex(reco_nu_vtx_sce_x,reco_nu_vtx_sce_y,reco_nu_vtx_sce_z);
	  bool in_fiducial_volume_reco = inFV(Reco_Vertex);


	  //topological score

	  if(topological_score > topo_cut){
	    sel_passed_topo_cut = true;
	  }


	  //end of topological score






	  int uncontained = 0;
	  int nonproton = 0;
	  int nPrimaryShowers = 0;
	  int nPrimaryTracks = 0;
	  std::vector<std::pair<int, float>> trk_llr_proto;
	  std::vector<std::pair<int, float>> trk_ilen;


	  for(int i = 0; i< trk_len_v->size(); i++) {
	    TVector3 endpoint(trk_sce_end_x_v->at(i),trk_sce_end_y_v->at(i),trk_sce_end_z_v->at(i));

	    if(pfp_generation_v->at(i) !=2) continue;

	    if(trk_score_v->at(i) >= 0.5){
	      nPrimaryTracks++;

	      if(!isContained(endpoint) && (trk_len_v->at(i) > 0) )  uncontained++;

	      if(trk_llr_pid_score_v->at(i) > 0.1  )  nonproton++;


	      trk_llr_proto.push_back({i, trk_llr_pid_score_v->at(i)});
	      trk_ilen.push_back({i, trk_len_v->at(i)});

	    }
	    else{
	      nPrimaryShowers++;
	      shower_index = i;
	    }

	  }



	  std::sort(trk_llr_proto.begin(), trk_llr_proto.end(), [](std::pair<int, float> a, std::pair<int, float> b) {
	    return a.second > b.second;
	  });


	  std::sort(trk_ilen.begin(), trk_ilen.end(), [](std::pair<int, float> a, std::pair<int, float> b) {
	    return a.second > b.second;
	  });



	  int i_max_len=0;
	  if(trk_ilen.size() > 0) i_max_len = trk_ilen.at(0).first; //index of the longest primary track



	  TVector3 reco_primary_vtx(reco_nu_vtx_sce_x,reco_nu_vtx_sce_y,reco_nu_vtx_sce_z);

	  TVector3 track_start(-1.0,-1.0,-1.0);
	  TVector3 nu_to_track_dist(-1.0,-1.0,-1.0);

	  if(trk_ilen.size() > 0){
	    track_start = TVector3 (trk_sce_start_x_v->at(i_max_len),trk_sce_start_y_v->at(i_max_len),trk_sce_start_z_v->at(i_max_len));
	    nu_to_track_dist = TVector3 (track_start - reco_primary_vtx);
	  }


	  for(size_t i_a = 0; i_a < trk_len_v->size(); i_a++){


	    if((pfp_generation_v->at(i_a) ==2) &&(trk_score_v->at(i_a) >= 0.5)
	       && (trk_llr_pid_score_v->at(i_a) > -1) && (trk_llr_pid_score_v->at(i_a) < 2)
	       && (trk_bragg_p_v->at(i_a)  > 0) && (trk_bragg_p_v->at(i_a)  < 500)
	       && (trk_bragg_mu_v->at(i_a) > 0) && (trk_bragg_mu_v->at(i_a) < 500)
	       && (trk_bragg_mip_v->at(i_a)> 0) && (trk_bragg_mip_v->at(i_a)< 500)
	       && (trk_len_v->at(i_a) < 1e6) && (trk_len_v->at(i_a) > 0)
	       ){

	      float tmvaOutput_mip = 0.0;
	      //	  float tmvaOutput_mu = 0.0;


	      trk_bragg_p_v_tmva = trk_bragg_p_v->at(i_a);
	      trk_bragg_mu_v_tmva = trk_bragg_mu_v->at(i_a);
	      trk_bragg_mip_v_tmva = trk_bragg_mip_v->at(i_a);
	      trk_llr_pid_score_v_tmva = trk_llr_pid_score_v->at(i_a);
	      trk_score_v_tmva = trk_score_v->at(i_a);
	      trk_len_v_tmva = trk_len_v->at(i_a);
	      trk_sce_end_x_v_tmva = trk_sce_end_x_v->at(i_a);
	      trk_sce_end_y_v_tmva = trk_sce_end_y_v->at(i_a);
	      trk_sce_end_z_v_tmva = trk_sce_end_z_v->at(i_a);
	      tmvaOutput_mip = tmvaReader->EvaluateMVA("BDT");



	      TVector3 track_start_ia (trk_sce_start_x_v->at(i_a),trk_sce_start_y_v->at(i_a),trk_sce_start_z_v->at(i_a));
	      TVector3 nu_to_track_dist_ia (track_start_ia - reco_primary_vtx);

	      //chaned dist to dist_ia to check if it makes ay difference.

	      if(  (trk_score_v->at(i_a) > muon_trk_score) && (nu_to_track_dist_ia.Mag() <= 4)
		   &&(trk_len_v->at(i_a) > muon_trk_len) && (trk_llr_pid_score_v->at(i_a) > muon_pid_score)
		   &&  (tmvaOutput_mip >= -0.1)){


		sel_has_muon_candidate = true;

		if(muon_index == -1){
		  muon_index = i_a;
		}
		else if ((trk_len_v->at(i_a) > trk_len_v->at(muon_index))){
		  muon_index = i_a;
		}

	      }
	    }

	  }

	  if(muon_index != -1){
	    muon_in_gap = ((pfnplanehits_U->at(muon_index)>=1)&&
			   (pfnplanehits_V->at(muon_index)>=1)&& (pfnplanehits_Y->at(muon_index)>=1));
	  }





	  //PION INDEX SELECTION

	  for (size_t i_b = 0; i_b < trk_len_v->size(); i_b++){

	    if((pfp_generation_v->at(i_b) ==2) &&(trk_score_v->at(i_b) >= 0.5)
	       && (trk_llr_pid_score_v->at(i_b) > -1) && (trk_llr_pid_score_v->at(i_b) < 2)
	       && (trk_bragg_p_v->at(i_b)  > 0) && (trk_bragg_p_v->at(i_b)  < 500)
	       && (trk_bragg_mu_v->at(i_b) > 0) && (trk_bragg_mu_v->at(i_b) < 500)
	       && (trk_bragg_mip_v->at(i_b)> 0) && (trk_bragg_mip_v->at(i_b)< 500)
	       && (trk_len_v->at(i_b) < 1e6) && (trk_len_v->at(i_b) > 0)
	       ){



	      TVector3 track_start_ib (trk_sce_start_x_v->at(i_b),trk_sce_start_y_v->at(i_b),trk_sce_start_z_v->at(i_b));
	      TVector3 nu_to_track_dist_ib (track_start_ib - reco_primary_vtx);




	      float tmvaOutput = 0.0;
	      float tmvaOutput_pi = 0.0;

	      trk_bragg_p_v_tmva = trk_bragg_p_v->at(i_b);
	      trk_bragg_mu_v_tmva = trk_bragg_mu_v->at(i_b);
	      trk_bragg_mip_v_tmva = trk_bragg_mip_v->at(i_b);
	      trk_llr_pid_score_v_tmva = trk_llr_pid_score_v->at(i_b);
	      trk_score_v_tmva = trk_score_v->at(i_b);
	      trk_len_v_tmva = trk_len_v->at(i_b);
	      trk_sce_end_x_v_tmva = trk_sce_end_x_v->at(i_b);
	      trk_sce_end_y_v_tmva = trk_sce_end_y_v->at(i_b);
	      trk_sce_end_z_v_tmva = trk_sce_end_z_v->at(i_b);
	      tmvaOutput = tmvaReader->EvaluateMVA("BDT");


	      trk_bragg_p_v_tmva_pi = trk_bragg_p_v->at(i_b);
	      trk_bragg_mu_v_tmva_pi = trk_bragg_mu_v->at(i_b);
	      trk_bragg_mip_v_tmva_pi = trk_bragg_mip_v->at(i_b);
	      trk_llr_pid_score_v_tmva_pi = trk_llr_pid_score_v->at(i_b);
	      trk_score_v_tmva_pi = trk_score_v->at(i_b);
	      trk_len_v_tmva_pi = trk_len_v->at(i_b);
	      trk_sce_end_x_v_tmva_pi = trk_sce_end_x_v->at(i_b);
	      trk_sce_end_y_v_tmva_pi = trk_sce_end_y_v->at(i_b);
	      trk_sce_end_z_v_tmva_pi = trk_sce_end_z_v->at(i_b);
	      tmvaOutput_pi = tmvaReader_pi->EvaluateMVA("BDT");



	      if ( (trk_llr_pid_score_v->at(i_b) > 0.1) && (nu_to_track_dist_ib.Mag() < 9.0)
		   && (i_b != muon_index)&& (tmvaOutput >= -0.1) && (tmvaOutput_pi >= -0.1)){


		pion_number++;
		pion_index = i_b;


		if(pion_index == -1){
		  pion_index = i_b;
		}
		else if ((trk_llr_pid_score_v->at(i_b) > trk_llr_pid_score_v->at(pion_index))){
		  pion_index = i_b;
		}
	      }
	    }
	  }


	  if(pion_index != -1){
	    pion_in_gap = ((pfnplanehits_U->at(pion_index)>0) && (pfnplanehits_V->at(pion_index)>0)&& (pfnplanehits_Y->at(pion_index)>0));
	    TVector3 pion_endpoint(trk_sce_end_x_v->at(pion_index),trk_sce_end_y_v->at(pion_index),trk_sce_end_z_v->at(pion_index));

	    if ((isContained(pion_endpoint)) )      	sel_contained_pions = true;

	  }




	  bool opening_angle_cut=false;
	  if((muon_index!= -1) && (pion_index!= -1)){
	    TVector3 MU(trk_dir_x_v->at(muon_index),trk_dir_y_v->at(muon_index),trk_dir_z_v->at(muon_index));
	    TVector3 PI(trk_dir_x_v->at(pion_index),trk_dir_y_v->at(pion_index),trk_dir_z_v->at(pion_index));
	    mu_pi_opening_angle = MU.Angle(PI);
	    if(mu_pi_opening_angle<2.6) opening_angle_cut = true;
	  }



	  //Shower cut


	  if(nPrimaryShowers==1){

	    TVector3 shower_start (trk_sce_start_x_v->at(shower_index),trk_sce_start_y_v->at(shower_index),trk_sce_start_z_v->at(shower_index));
	    TVector3 nu_to_shower_dist (shower_start - reco_primary_vtx);



	    // std::cout<<nPrimaryShowers<<"\t"<<shower_index<<"\t"<<nu_to_shower_dist.Mag()<<std::endl;

	    if((pfnplanehits_Y->at(shower_index)<50) && (pfnplanehits_Y->at(shower_index)>test) && (nu_to_shower_dist.Mag()>10)){

	      shower_cut = true ;
	    }

	    //	  if(shr_moliere_avg_v->at(i_c) < 200){

	    //	    mol_avg = true;
	    //	  }

	  }
	  else if (nPrimaryShowers==0)  shower_cut = true ;



          int rp = pion_number;// number of reconstructed pions for the selection cuts


	  Cuts[0][rp] = true;// in_fiducial_volume_true;  //all events with entries in true FV
	  Cuts[1][rp] = Cuts[0][rp] && in_fiducial_volume_reco; // Events in fiducial volume
	  Cuts[2][rp] = Cuts[1][rp] && sel_passed_topo_cut;// Topological score
	  Cuts[3][rp] = Cuts[2][rp] && sel_has_muon_candidate;// Muon candidate
	  Cuts[4][rp] = Cuts[3][rp] && sel_contained_pions; //Contained charged pion
	  Cuts[5][rp] = Cuts[4][rp] && muon_in_gap;
	  Cuts[6][rp] = Cuts[5][rp] && pion_in_gap;
	  Cuts[7][rp] = Cuts[6][rp] && shower_cut; //no or one shower with conditions
	  Cuts[8][rp] = Cuts[7][rp] && opening_angle_cut; //muon-pion opening angle cut
	  Cuts[9][rp] = Cuts[8][rp] && (nonproton == 2); // no other charged pions





	  //Filling histograms for all variables and cuts for signal and background
	  for(int c=0; c<ncuts; c++){
	     for(int p=0;p<npions;p++){
	    int s=-1;
	    if((i_f < 4) && signal[p] ) s = 0;
	    else if ((i_f <4) && !signal[p]) s = 1;
	    else if(i_f == 4) s = 2;//EXT
	    else if(i_f == 5) s = 3;//Dirt
	    else if(i_f == 6) s = 4;//Data
	    else continue;

	    if (Cuts[c][p]){

	      Selected[c][s][p] += Scale[i_f];
//	      std::cout<<CutsName[c]<<"\t"<<Sample[s]<<"\t"<<i_f<<"\t"<<  Scale[i_f]<<"\t"<<Selected[c][s]<<std::endl;

	      Histos[c][s][0][p]->Fill(topological_score,Scale[i_f]); //Variable[0]="Topological Score";

	      if(muon_index!=-1) {
		TVector3 muon_start (trk_sce_start_x_v->at(muon_index),trk_sce_start_y_v->at(muon_index),trk_sce_start_z_v->at(muon_index));
		TVector3 nu_to_track_dist (muon_start - reco_primary_vtx);

		Histos[c][s][1][p]->Fill(trk_len_v->at(muon_index),Scale[i_f] ); //Variable[1]="Muon track length";
		Histos[c][s][4][p]->Fill(trk_llr_pid_score_v->at(muon_index), Scale[i_f] );//  Variable[4]="Muon PID";
		Histos[c][s][6][p]->Fill(nu_to_track_dist.Mag(), Scale[i_f]); //Variable[6]="MuonVtxDistance";
		Histos[c][s][9][p]->Fill(pfnplanehits_U->at(muon_index), Scale[i_f]);  //Variable[9]="MuonUPlaneHits";
		Histos[c][s][10][p]->Fill(pfnplanehits_Y->at(muon_index), Scale[i_f]); //Variable[10]="MuonYPlaneHits";
		Histos[c][s][11][p]->Fill(pfnplanehits_V->at(muon_index), Scale[i_f]); //Variable[11]="MuonVPlaneHits";

		float tmvaOutput_mip = 0.0;

		int i_a=muon_index;
		trk_bragg_p_v_tmva = trk_bragg_p_v->at(i_a);
		trk_bragg_mu_v_tmva = trk_bragg_mu_v->at(i_a);
		trk_bragg_mip_v_tmva = trk_bragg_mip_v->at(i_a);
		trk_llr_pid_score_v_tmva = trk_llr_pid_score_v->at(i_a);
		trk_score_v_tmva = trk_score_v->at(i_a);
		trk_len_v_tmva = trk_len_v->at(i_a);
		trk_sce_end_x_v_tmva = trk_sce_end_x_v->at(i_a);
		trk_sce_end_y_v_tmva = trk_sce_end_y_v->at(i_a);
		trk_sce_end_z_v_tmva = trk_sce_end_z_v->at(i_a);
		tmvaOutput_mip = tmvaReader->EvaluateMVA("BDT");

		Histos[c][s][23][p]->Fill(tmvaOutput_mip, Scale[i_f] );  //Variable[23]="PionTMVAMip";


	      }

	      if(pion_index!=-1) {
		TVector3 pion_start (trk_sce_start_x_v->at(pion_index),trk_sce_start_y_v->at(pion_index),trk_sce_start_z_v->at(pion_index));
		TVector3 nu_to_track_dist (pion_start - reco_primary_vtx);

		Histos[c][s][2][p]->Fill(trk_len_v->at(pion_index), Scale[i_f]);//  Variable[2]="Pion track length";
		Histos[c][s][5][p]->Fill(trk_llr_pid_score_v->at(pion_index),Scale[i_f] );//  Variable[5]="Pion PID";
		Histos[c][s][7][p]->Fill(nu_to_track_dist.Mag(),Scale[i_f]); //Variable[7]="PionVtxDistance";
		Histos[c][s][12][p]->Fill(pfnplanehits_U->at(pion_index),Scale[i_f]);  //Variable[12]="PionUPlaneHits";
		Histos[c][s][13][p]->Fill(pfnplanehits_Y->at(pion_index),Scale[i_f]);  //Variable[13]="PionYPlaneHits";
		Histos[c][s][14][p]->Fill(pfnplanehits_V->at(pion_index),Scale[i_f]);  //Variable[14]="PionVPlaneHits";
		float tmvaOutput = 0.0;
		float tmvaOutput_pi = 0.0;
		int i_b=pion_index;
		trk_bragg_p_v_tmva = trk_bragg_p_v->at(i_b);
		trk_bragg_mu_v_tmva = trk_bragg_mu_v->at(i_b);
		trk_bragg_mip_v_tmva = trk_bragg_mip_v->at(i_b);
		trk_llr_pid_score_v_tmva = trk_llr_pid_score_v->at(i_b);
		trk_score_v_tmva = trk_score_v->at(i_b);
		trk_len_v_tmva = trk_len_v->at(i_b);
		trk_sce_end_x_v_tmva = trk_sce_end_x_v->at(i_b);
		trk_sce_end_y_v_tmva = trk_sce_end_y_v->at(i_b);
		trk_sce_end_z_v_tmva = trk_sce_end_z_v->at(i_b);
		tmvaOutput = tmvaReader->EvaluateMVA("BDT");


		trk_bragg_p_v_tmva_pi = trk_bragg_p_v->at(i_b);
		trk_bragg_mu_v_tmva_pi = trk_bragg_mu_v->at(i_b);
		trk_bragg_mip_v_tmva_pi = trk_bragg_mip_v->at(i_b);
		trk_llr_pid_score_v_tmva_pi = trk_llr_pid_score_v->at(i_b);
		trk_score_v_tmva_pi = trk_score_v->at(i_b);
		trk_len_v_tmva_pi = trk_len_v->at(i_b);
		trk_sce_end_x_v_tmva_pi = trk_sce_end_x_v->at(i_b);
		trk_sce_end_y_v_tmva_pi = trk_sce_end_y_v->at(i_b);
		trk_sce_end_z_v_tmva_pi = trk_sce_end_z_v->at(i_b);
		tmvaOutput_pi = tmvaReader_pi->EvaluateMVA("BDT");


		Histos[c][s][21][p]->Fill(tmvaOutput,Scale[i_f] );  //Variable[21]="PionTMVAMip";
		Histos[c][s][22][p]->Fill(tmvaOutput_pi,Scale[i_f] );  //Variable[22]="PionTMVAPi";


	      }
	      if((muon_index!= -1) && (pion_index!= -1)) {
		TVector3 muon_start (trk_sce_start_x_v->at(muon_index),trk_sce_start_y_v->at(muon_index),trk_sce_start_z_v->at(muon_index));
		TVector3 pion_start (trk_sce_start_x_v->at(pion_index),trk_sce_start_y_v->at(pion_index),trk_sce_start_z_v->at(pion_index));
		TVector3 pion_muon_dist (muon_start - pion_start);


		Histos[c][s][3][p]->Fill(mu_pi_opening_angle,Scale[i_f]);//  Variable[3]="Mu-Pi Opening angle";
		Histos[c][s][8][p]->Fill(pion_muon_dist.Mag(),Scale[i_f]); // Variable[8]="MuonPionDistance";














	      }
	      if(shower_index !=-1){

		TVector3 shower_start (trk_sce_start_x_v->at(shower_index),trk_sce_start_y_v->at(shower_index),trk_sce_start_z_v->at(shower_index));
		TVector3 nu_to_shower_dist (shower_start - reco_primary_vtx);


		Histos[c][s][15][p]->Fill(pfnplanehits_U->at(shower_index),Scale[i_f]);  //Variable[15]="ShowerUPlaneHits";
		Histos[c][s][16][p]->Fill(pfnplanehits_Y->at(shower_index),Scale[i_f]);  //Variable[16]="ShowerYPlaneHits";
		Histos[c][s][17][p]->Fill(pfnplanehits_V->at(shower_index),Scale[i_f]);  //Variable[17]="ShowerVPlaneHits";

		Histos[c][s][20][p]->Fill(nu_to_shower_dist.Mag(),Scale[i_f]);// Variable[20]="ShowerVtxDistance";

	      }

	      Histos[c][s][18][p]->Fill(nPrimaryShowers,Scale[i_f]); // Variable[18]="NPrimaryShowers";
	      Histos[c][s][19][p]->Fill(nPrimaryTracks,Scale[i_f]); // Variable[19]="NPrimaryTracks";





	    }
	     }
	  }





	} //event loop

      } //for files

      output.close();


      for(int p=0;p<npions;p++){

      TFile *myfile = new TFile("histograms.root","RECREATE");
      for(int c=0;c<ncuts;c++){
        for(int v=0;v<nvariables;v++){

	  for(int s=0;s<nsamples;s++){

	    Histos[c][s][v][p]->Write();
	  }
        }
      }




	
      double purity = 0.0;
      double efficiency = 0.0;
      std::cout<<"pions="<<p<<"\t \t \t \t signal  \t \t bckg \t \t \t EXT \t \t \t Dirt \t \t \t DATA \t \t \t efficiency \t purity \t eff*pur\n";
      for(int c=0; c<ncuts;c++){

	double purity = Selected[c][0][p]/double(Selected[c][0][p]+Selected[c][1][p] + Selected[c][2][p] +Selected[c][3][p] )*100 ;
	double efficiency = Selected[c][0][p]/double(Selected[0][0][p])*100;

	std::cout<<CutsName[c]<<"\t \t \t "<<Selected[c][0][p]
		 <<"\t \t "<<Selected[c][1][p]
		 <<"\t \t "<<Selected[c][2][p]
		 <<"\t \t "<<Selected[c][3][p]
		 <<"\t \t "<<Selected[c][4][p]
		 <<"\t\t\t"<<efficiency
		 <<"\t \t "<<purity
		 <<"\t\t\t"<<efficiency*purity
		 << std::endl;

      }


      /*
      TCanvas *c1 = new TCanvas("c1","",2000,800);
      c1->Divide(5,2);

      for(int i = 0; i<nvariables; i++){
	for(int c=0;c<ncuts;c++){

	  c1->cd(c+1);



	  Histos[c][0][i][p]->Draw();
	  Histos[c][1][i][p]->Draw("same");
	  Histos[c][2][i][p]->Draw("same");
	  Histos[c][3][i][p]->Draw("same");
	  Histos[c][4][i][p]->Draw("same e0");
	  // Histos[c][0][i][p]->Draw("same");


	}
	c1->Print("plots/"+Variable[i]+".png");
      }
*/


      }//pions 

    }//TEST END
}
