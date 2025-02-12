The existing script is a procedural implementation for analyzing neutrino interaction events. We'll refactor it into an object-oriented form by creating a class that encapsulates the functionality and the data members.

Here's a possible object-oriented version of the script:

```c++
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
#include <vector>
#include <string>
#include <fstream>

class CCpiAnalysis {
public:
    CCpiAnalysis();
    void runAnalysis();

private:
    const int ncuts = 10;
    const int nvariables = 26;
    const int nsamples = 6;

    double Selected[10][6] = {0};
    TString CutsName[10];
    TString Variable[26];
    TString Sample[6];
    TH1F *Histos[10][6][26];
    std::vector<std::string> Files;
    double dataPoT[4] = {0};
    double MCPoT[4] = {0};
    double Scale[7] = {1};

    TMVA::Reader *tmvaReader;
    TMVA::Reader *tmvaReader_mu;
    TMVA::Reader *tmvaReader_pi;

    void initialize();
    void processFile(const std::string& filename, int fileIndex);
    bool inFV(const TVector3& vertex);
    bool isContained(const TVector3& endpoint);
};

CCpiAnalysis::CCpiAnalysis() {
    initialize();
}

void CCpiAnalysis::initialize() {
    CutsName[0] = "EventsInTrueFV";
    CutsName[1] = "InFiducialVol";
    CutsName[2] = "Topological";
    CutsName[3] = "MuonCandidate";
    CutsName[4] = "ContainedPion";
    CutsName[5] = "MuonIn3Planes";
    CutsName[6] = "PionIn3Planes";
    CutsName[7] = "ShowerCut";
    CutsName[8] = "OpeningAngle";
    CutsName[9] = "Nonprotons";

    Variable[0] = "TopologicalScore";
    Variable[1] = "MuonTrackLength";
    Variable[2] = "PionTrackLength";
    Variable[3] = "MuPiOpeningangle";
    Variable[4] = "MuonPID";
    Variable[5] = "PionPID";
    Variable[6] = "MuonVtxDistance";
    Variable[7] = "PionVtxDistance";
    Variable[8] = "MuonPionDistance";
    Variable[9] = "MuonUPlaneHits";
    Variable[10] = "MuonYPlaneHits";
    Variable[11] = "MuonVPlaneHits";
    Variable[12] = "PionUPlaneHits";
    Variable[13] = "PionYPlaneHits";
    Variable[14] = "PionVPlaneHits";
    Variable[15] = "ShowerUPlaneHits";
    Variable[16] = "ShowerYPlaneHits";
    Variable[17] = "ShowerVPlaneHits";
    Variable[18] = "NPrimaryShowers";
    Variable[19] = "NPrimaryTracks";
    Variable[20] = "ShowerVtxDistance";
    Variable[21] = "PionTMVAMip";
    Variable[22] = "PionTMVAPi";
    Variable[23] = "MuonTMVAMip";
    Variable[24] = "PionTrkEnergy";
    Variable[25] = "MuonTrkEnrgy";

    Sample[0] = "Signal";
    Sample[1] = "Background";
    Sample[2] = "EXT";
    Sample[3] = "OOFV";
    Sample[4] = "Data";
    Sample[5] = "MC";

    for (int c = 0; c < ncuts; ++c) {
        for (int v = 0; v < nvariables; ++v) {
            for (int s = 0; s < nsamples; ++s) {
                Histos[c][s][v] = new TH1F(Variable[v] + "_" + CutsName[c] + "_" + Sample[s], "", 100, -1, -1);
                if (s != 4) Histos[c][s][v]->SetLineColor(s + 1);
                Histos[c][s][v]->GetXaxis()->SetTitle(Variable[v]);
                Histos[c][s][v]->SetTitle(CutsName[c]);
            }
        }
    }

    Files = {
        "/data/uboone/new_numi_flux/Run1_fhc_new_numi_flux_fhc_pandora_ntuple.root",
        "/data/uboone/new_numi_flux/Run2_fhc_new_numi_flux_fhc_pandora_ntuple.root",
        "/data/uboone/new_numi_flux/Run4_fhc_new_numi_flux_fhc_pandora_ntuple.root",
        "/data/uboone/new_numi_flux/Run5_fhc_new_numi_flux_fhc_pandora_ntuple.root",
        "/data/uboone/EXT/neutrinoselection_filt_run1_beamoff.root",
        "/data/uboone/dirt/prodgenie_numi_uboone_overlay_dirt_fhc_mcc9_run1_v28_all_snapshot.root"
    };

    dataPoT[0] = 2.192 + 3.283; // Run1 including open trigger
    dataPoT[1] = 1.268;  // Run2
    dataPoT[2] = 2.075;  // Run4
    dataPoT[3] = 2.231;  // Run5

    MCPoT[0] = 23.282;  // Run1
    MCPoT[1] = 24.9337; // Run2
    MCPoT[2] = 11.04 + 17.295; // Run4
    MCPoT[3] = 63.145; // Run5

    double totalDataPoT = 0;
    double totalMCPoT = 0;

    for (int i = 0; i < 4; ++i) {
        Scale[i] = dataPoT[i] / MCPoT[i];      // data/simulation;
        totalDataPoT += dataPoT[i];
        totalMCPoT += MCPoT[i];
    }

    Scale[5] = totalDataPoT / 16.739;      // Run1 dirt from David 1.67392e+21
    Scale[4] = 5748692.0 / 9199232.0;  // Run1 beam off triggers scaling to beam on triggers

    tmvaReader = new TMVA::Reader();
    tmvaReader_mu = new TMVA::Reader();
    tmvaReader_pi = new TMVA::Reader();

    tmvaReader->AddVariable("trk_bragg_p_v", &trk_bragg_p_v_tmva);
    tmvaReader->AddVariable("trk_bragg_mu_v", &trk_bragg_mu_v_tmva);
    tmvaReader->AddVariable("trk_bragg_mip_v", &trk_bragg_mip_v_tmva);
    tmvaReader->AddVariable("trk_llr_pid_score_v", &trk_llr_pid_score_v_tmva);
    tmvaReader->AddVariable("trk_score_v", &trk_score_v_tmva);
    tmvaReader->AddVariable("trk_len_v", &trk_len_v_tmva);
    tmvaReader->AddVariable("trk_sce_end_x_v", &trk_sce_end_x_v_tmva);
    tmvaReader->AddVariable("trk_sce_end_y_v", &trk_sce_end_y_v_tmva);
    tmvaReader->AddVariable("trk_sce_end_z_v", &trk_sce_end_z_v_tmva);

    tmvaReader->BookMVA("BDT", "booster_decision_tree/dataset_MIP_BDT/weights/TMVAClassification_BDT.weights.xml");

    tmvaReader_mu->AddVariable("trk_bragg_p_v", &trk_bragg_p_v_tmva_mu);
    tmvaReader_mu->AddVariable("trk_bragg_mu_v", &trk_bragg_mu_v_tmva_mu);
    tmvaReader_mu->AddVariable("trk_bragg_mip_v", &trk_bragg_mip_v_tmva_mu);
    tmvaReader_mu->AddVariable("trk_llr_pid_score_v", &trk_llr_pid_score_v_tmva_mu);
    tmvaReader_mu->AddVariable("trk_score_v", &trk_score_v_tmva_mu);
    tmvaReader_mu->AddVariable("trk_len_v", &trk_len_v_tmva_mu);
    tmvaReader_mu->AddVariable("trk_sce_end_x_v", &trk_sce_end_x_v_tmva_mu);
    tmvaReader_mu->AddVariable("trk_sce_end_y_v", &trk_sce_end_y_v_tmva_mu);
    tmvaReader_mu->AddVariable("trk_sce_end_z_v", &trk_sce_end_z_v_tmva_mu);

    tmvaReader_mu->BookMVA("BDT", "booster_decision_tree/dataset_muon_BDT/weights/TMVAClassification_BDT.weights.xml");

    tmvaReader_pi->AddVariable("trk_bragg_p_v", &trk_bragg_p_v_tmva_pi);
    tmvaReader_pi->AddVariable("trk_bragg_mu_v", &trk_bragg_mu_v_tmva_pi);
    tmvaReader_pi->AddVariable("trk_bragg_mip_v", &trk_bragg_mip_v_tmva_pi);
    tmvaReader_pi->AddVariable("trk_llr_pid_score_v", &trk_llr_pid_score_v_tmva_pi);
    tmvaReader_pi->AddVariable("trk_score_v", &trk_score_v_tmva_pi);
    tmvaReader_pi->AddVariable("trk_len_v", &trk_len_v_tmva_pi);
    tmvaReader_pi->AddVariable("trk_sce_end_x_v", &trk_sce_end_x_v_tmva_pi);
    tmvaReader_pi->AddVariable("trk_sce_end_y_v", &trk_sce_end_y_v_tmva_pi);
    tmvaReader_pi->AddVariable("trk_sce_end_z_v", &trk_sce_end_z_v_tmva_pi);

    tmvaReader_pi->BookMVA("BDT", "booster_decision_tree/dataset_pion_BDT/weights/TMVAClassification_BDT.weights.xml");
}

void CCpiAnalysis::runAnalysis() {
    std::ofstream output("backtracked_showers.txt");

    double test_min = 0.0;
    double test_max = 0.03;
    double test_step = 0.1;

    for (double test = test_min; test < test_max; test += test_step) {
        for (int cc = 0; cc < ncuts; cc++) {
            for (int j = 0; j < 7; j++) {
                Selected[cc][j] = 0;
            }
        }

        for (size_t i_f = 0; i_f < Files.size(); i_f++) {
            processFile(Files[i_f], i_f);
        }

        output.close();

        TFile *myfile = new TFile("histograms.root", "RECREATE");
        for (int c = 0; c < ncuts; c++) {
            for (int v = 0; v < nvariables; v++) {
                for (int s = 0; s < nsamples; s++) {
                    Histos[c][s][v]->Write();
                }
            }
        }

        double purity = 0.0;
        double efficiency = 0.0;
        std::cout << "test=" << test << "\t \t \t \t signal  \t \t bckg \t \t \t EXT \t \t \t Dirt \t \t \t DATA \t \t \t efficiency \t purity \t eff*pur\n";
        for (int c = 0; c < ncuts; c++) {
            purity = Selected[c][0] / double(Selected[c][0] + Selected[c][1] + Selected[c][2] + Selected[c][3]) * 100;
            efficiency = Selected[c][0] / double(Selected[0][0]) * 100;

            std::cout << CutsName[c] << "\t \t \t " << Selected[c][0] << "(" << Selected[c][0] << ")"
                      << "\t \t " << Selected[c][1] << "(" << Selected[c][1] << ")"
                      << "\t \t " << Selected[c][2] << "(" << Selected[c][2] << ")"
                      << "\t \t " << Selected[c][3] << "(" << Selected[c][3] << ")"
                      << "\t \t " << Selected[c][4]
                      << "\t\t\t" << efficiency
                      << "\t \t " << purity
                      << "\t\t\t" << efficiency * purity
                      << std::endl;
        }

        TCanvas *c1 = new TCanvas("c1", "", 2000, 800);
        c1->Divide(5, 2);

        for (int i = 0; i < nvariables; i++) {
            for (int c = 0; c < ncuts; c++) {
                c1->cd(c + 1);

                Histos[c][0][i]->Draw();
                Histos[c][1][i]->Draw("same");
                Histos[c][2][i]->Draw("same");
                Histos[c][3][i]->Draw("same");
                Histos[c][4][i]->Draw("same e0");
            }
            c1->Print("plots/" + Variable[i] + ".png");
        }
    }
}

bool CCpiAnalysis::inFV(const TVector3& vertex) {
    // Implement the inFV function here
    return true; // Placeholder implementation
}

bool CCpiAnalysis::isContained(const TVector3& endpoint) {
    // Implement the isContained function here
    return true; // Placeholder implementation
}

void CCpiAnalysis::processFile(const std::string& filename, int fileIndex) {
    TFile* f = TFile::Open(filename.c_str());
    if ((!f) || f->IsZombie()) {
        delete f;
        return;
    }
    TTree *t;
    f->GetObject("nuselection/NeutrinoSelectionFilter", t);

    if (!t) {
        delete f;
        std::cout << "No tree NeutrinoSelectionFilter found in file " << filename << std::endl;
        return;
    }

    t->SetMakeClass(1);
    t->SetBranchStatus("*", 0);

    Int_t run, sub, evt, nu_pdg, ccnc, interaction;
    Float_t true_nu_px, true_nu_py, true_nu_pz, reco_nu_vtx_sce_x, reco_nu_vtx_sce_y, reco_nu_vtx_sce_z;
    std::vector<int> *backtracked_pdg = nullptr;
    std::vector<float> *backtracked_start_x = nullptr;
    std::vector<float> *backtracked_start_y = nullptr;
    std::vector<float> *backtracked_start_z = nullptr;
    Int_t n_tracks, n_showers;
    std::vector<float> *trk_score_v = nullptr;
    std::vector<int> *mc_pdg = nullptr;
    std::vector<float> *mc_E = nullptr;
    std::vector<float> *mc_vx = nullptr;
    std::vector<float> *mc_vy = nullptr;
    std::vector<float> *mc_vz = nullptr;
    std::vector<float> *mc_endx = nullptr;
    std::vector<float> *mc_endy = nullptr;
    std::vector<float> *mc_endz = nullptr;
    std::vector<float> *mc_px = nullptr;
    std::vector<float> *mc_py = nullptr;
    std::vector<float> *mc_pz = nullptr;
    std::vector<float> *trk_dir_x_v = nullptr;
    std::vector<float> *trk_dir_y_v = nullptr;
    std::vector<float> *trk_dir_z_v = nullptr;
    std::vector<float> *trk_sce_start_x_v = nullptr;
    std::vector<float> *trk_sce_start_y_v = nullptr;
    std::vector<float> *trk_sce_start_z_v = nullptr;
    std::vector<float> *trk_sce_end_x_v = nullptr;
    std::vector<float> *trk_sce_end_y_v = nullptr;
    std::vector<float> *trk_sce_end_z_v = nullptr;
    std::vector<float> *trk_len_v = nullptr;
    std::vector<float> *trk_calo_energy_u_v = nullptr;
    std::vector<float> *trk_calo_energy_v_v = nullptr;
    std::vector<float> *trk_calo_energy_y_v = nullptr;
    std::vector<float> *trk_llr_pid_score_v = nullptr;
    Float_t topological_score;
    std::vector<int> *pfnplanehits_U = nullptr;
    std::vector<int> *pfnplanehits_V = nullptr;
    std::vector<int> *pfnplanehits_Y = nullptr;
    std::vector<unsigned int> *pfp_generation_v = nullptr;
    std::vector<float> *trk_bragg_p_v = nullptr;
    std::vector<float> *trk_bragg_mu_v = nullptr;
    std::vector<float> *trk_bragg_mip_v = nullptr;
    std::vector<float> *shr_start_x_v = nullptr;
    std::vector<float> *shr_start_y_v = nullptr;
    std::vector<float> *shr_start_z_v = nullptr;
    std::vector<float> *shr_dist_v = nullptr;
    std::vector<float> *shr_moliere_avg_v = nullptr;
    std::vector<float> *shr_moliere_rms_v = nullptr;
    std::vector<float> *shr_tkfit_dedx_u_v = nullptr;
    std::vector<float> *shr_tkfit_dedx_v_v = nullptr;
    std::vector<float> *shr_tkfit_dedx_y_v = nullptr;

    t->SetBranchAddress("run", &run);
    t->SetBranchAddress("sub", &sub);
    t->SetBranchAddress("evt", &evt);
    t->SetBranchAddress("nu_pdg", &nu_pdg);
    t->SetBranchAddress("ccnc", &ccnc);
    t->SetBranchAddress("interaction", &interaction);
    t->SetBranchAddress("true_nu_px", &true_nu_px);
