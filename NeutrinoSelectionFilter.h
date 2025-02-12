//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Feb 12 11:33:38 2025 by ROOT version 6.34.02
// from TTree NeutrinoSelectionFilter/Neutrino Selection TTree
// found on file: Run1_fhc_new_numi_flux_fhc_pandora_ntuple.root
//////////////////////////////////////////////////////////

#ifndef NeutrinoSelectionFilter_h
#define NeutrinoSelectionFilter_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"
#include "vector"
#include "vector"
#include "string"
#include "vector"
#include "vector"
#include "vector"

class NeutrinoSelectionFilter {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           selected;
   Int_t           run;
   Int_t           sub;
   Int_t           evt;
   UInt_t          trk_id;
   UInt_t          shr_id;
   UInt_t          trk2_id;
   UInt_t          shr2_id;
   UInt_t          trk3_id;
   UInt_t          shr3_id;
   Float_t         shr_energy_tot;
   Float_t         shr_energy;
   Float_t         shr_energy_second;
   Float_t         shr_energy_third;
   Float_t         shr_energy_tot_cali;
   Float_t         shr_energy_cali;
   Float_t         shr_energy_second_cali;
   Float_t         shr_energy_third_cali;
   Float_t         shr_theta;
   Float_t         shr_phi;
   Float_t         shr_pca_0;
   Float_t         shr_pca_1;
   Float_t         shr_pca_2;
   Float_t         shr_px;
   Float_t         shr_py;
   Float_t         shr_pz;
   Float_t         shr_openangle;
   Float_t         shr_tkfit_start_x;
   Float_t         shr_tkfit_start_y;
   Float_t         shr_tkfit_start_z;
   Float_t         shr_tkfit_theta;
   Float_t         shr_tkfit_phi;
   Float_t         shr_start_x;
   Float_t         shr_start_y;
   Float_t         shr_start_z;
   Float_t         shr_dedx_Y;
   Float_t         shr_dedx_V;
   Float_t         shr_dedx_U;
   Float_t         shr_dedx_Y_cali;
   Float_t         shr_dedx_V_cali;
   Float_t         shr_dedx_U_cali;
   Float_t         shr_tkfit_dedx_Y;
   Float_t         shr_tkfit_dedx_V;
   Float_t         shr_tkfit_dedx_U;
   Float_t         shr_tkfit_dedx_max;
   UInt_t          shr_tkfit_nhits_Y;
   UInt_t          shr_tkfit_nhits_V;
   UInt_t          shr_tkfit_nhits_U;
   Float_t         shr_llrpid_dedx_Y;
   Float_t         shr_llrpid_dedx_V;
   Float_t         shr_llrpid_dedx_U;
   Float_t         shr_llrpid_dedx;
   Float_t         shr_tkfit_dedx_Y_alt;
   Float_t         shr_tkfit_dedx_V_alt;
   Float_t         shr_tkfit_dedx_U_alt;
   UInt_t          shr_tkfit_nhits_Y_alt;
   UInt_t          shr_tkfit_nhits_V_alt;
   UInt_t          shr_tkfit_nhits_U_alt;
   Float_t         trkfit;
   UInt_t          shr_tkfit_npoints;
   UInt_t          shr_tkfit_npointsvalid;
   Float_t         shr_trkfitmedangle;
   Float_t         shrmoliereavg;
   Float_t         shrmoliererms;
   Float_t         shr1shr2moliereavg;
   Float_t         shr1shr2moliererms;
   Float_t         shr1trk1moliereavg;
   Float_t         shr1trk1moliererms;
   Float_t         shr1trk2moliereavg;
   Float_t         shr1trk2moliererms;
   UChar_t         ismerged;
   Float_t         merge_bestdot;
   Float_t         merge_bestdist;
   Float_t         merge_vtx_x;
   Float_t         merge_vtx_y;
   Float_t         merge_vtx_z;
   UInt_t          merge_tk_ipfp;
   Float_t         shr_tkfit_2cm_dedx_Y;
   Float_t         shr_tkfit_2cm_dedx_V;
   Float_t         shr_tkfit_2cm_dedx_U;
   UInt_t          shr_tkfit_2cm_nhits_Y;
   UInt_t          shr_tkfit_2cm_nhits_V;
   UInt_t          shr_tkfit_2cm_nhits_U;
   Float_t         shr_tkfit_gap05_dedx_Y;
   Float_t         shr_tkfit_gap05_dedx_V;
   Float_t         shr_tkfit_gap05_dedx_U;
   UInt_t          shr_tkfit_gap05_nhits_Y;
   UInt_t          shr_tkfit_gap05_nhits_V;
   UInt_t          shr_tkfit_gap05_nhits_U;
   Float_t         shr_tkfit_gap10_dedx_Y;
   Float_t         shr_tkfit_gap10_dedx_V;
   Float_t         shr_tkfit_gap10_dedx_U;
   UInt_t          shr_tkfit_gap10_nhits_Y;
   UInt_t          shr_tkfit_gap10_nhits_V;
   UInt_t          shr_tkfit_gap10_nhits_U;
   Float_t         shr_chipr;
   Float_t         shr_chimu;
   Float_t         shr_bragg_p;
   Float_t         shr_bragg_mu;
   Float_t         shr_bragg_mip;
   Float_t         shr_bragg_kaon;
   Float_t         shr_bragg_pion;
   Float_t         tksh_distance;
   Float_t         tksh_angle;
   Float_t         shr_distance;
   Float_t         shr_score;
   Int_t           shr_bkt_pdg;
   Float_t         shr_bkt_purity;
   Float_t         shr_bkt_completeness;
   Float_t         shr_bkt_E;
   Float_t         trk_len;
   Float_t         trk_theta;
   Float_t         trk_phi;
   Float_t         trk_energy;
   Float_t         trk_energy_muon;
   Float_t         trk_energy_muon_mcs;
   Float_t         trk_energy_tot;
   Float_t         trk_energy_muon_tot;
   Float_t         trk_distance;
   Float_t         trk_score;
   Int_t           trk_bkt_pdg;
   Float_t         trk_bkt_purity;
   Float_t         trk_bkt_completeness;
   Float_t         trk_bkt_E;
   Float_t         trk_chipr_best;
   Float_t         trk_chipr_worst;
   Float_t         trk_chimu_best;
   Float_t         trk_chimu_worst;
   Float_t         trk_chipr;
   Float_t         trk_chimu;
   Float_t         trk_pida;
   Float_t         trk_bragg_p;
   Float_t         trk_bragg_mu;
   Float_t         trk_bragg_mip;
   Float_t         trk_bragg_kaon;
   Float_t         trk_bragg_pion;
   UInt_t          trk_hits_max;
   UInt_t          shr_hits_max;
   vector<int>     *all_shr_hits;
   vector<int>     *all_trk_hits;
   vector<float>   *all_shr_energies;
   vector<float>   *all_trk_energies;
   UInt_t          shr_hits_max;
   UInt_t          trk_hits_2nd;
   UInt_t          shr_hits_2nd;
   UInt_t          trk_hits_3rd;
   UInt_t          shr_hits_3rd;
   Float_t         trkshrhitdist0;
   Float_t         trkshrhitdist1;
   Float_t         trkshrhitdist2;
   Float_t         trk2shrhitdist0;
   Float_t         trk2shrhitdist1;
   Float_t         trk2shrhitdist2;
   Float_t         trk1trk2hitdist0;
   Float_t         trk1trk2hitdist1;
   Float_t         trk1trk2hitdist2;
   UInt_t          total_hits_y;
   Float_t         extra_energy_y;
   Float_t         trk_energy_hits_tot;
   UInt_t          subcluster;
   UInt_t          shrsubclusters0;
   UInt_t          shrsubclusters1;
   UInt_t          shrsubclusters2;
   Float_t         shrclusfrac0;
   Float_t         shrclusfrac1;
   Float_t         shrclusfrac2;
   Float_t         shrclusdir0;
   Float_t         shrclusdir1;
   Float_t         shrclusdir2;
   UInt_t          shr_hits_tot;
   UInt_t          shr_hits_y_tot;
   UInt_t          shr_hits_u_tot;
   UInt_t          shr_hits_v_tot;
   UInt_t          trk_hits_tot;
   UInt_t          trk_hits_y_tot;
   UInt_t          trk_hits_u_tot;
   UInt_t          trk_hits_v_tot;
   Float_t         _elecclusters_U_charge;
   Float_t         _elecclusters_V_charge;
   Float_t         _elecclusters_Y_charge;
   Int_t           _elecclusters_U_N;
   Int_t           _elecclusters_V_N;
   Int_t           _elecclusters_Y_N;
   UInt_t          n_tracks_contained;
   UInt_t          n_showers_contained;
   Float_t         matched_E;
   Float_t         hits_ratio;
   Float_t         contained_fraction;
   Float_t         sps_contained_fraction;
   Float_t         pt;
   Float_t         p;
   Float_t         pt_assume_muon;
   Float_t         p_assume_muon;
   Float_t         reco_e;
   Float_t         dvtx;
   Float_t         dtrk;
   Float_t         contained_sps_ratio;
   vector<double>  *dtrk_x_boundary;
   vector<double>  *dtrk_y_boundary;
   vector<double>  *dtrk_z_boundary;
   vector<double>  *dshr_x_boundary;
   vector<double>  *dshr_y_boundary;
   vector<double>  *dshr_z_boundary;
   vector<double>  *dvtx_x_boundary;
   vector<double>  *dvtx_y_boundary;
   vector<double>  *dvtx_z_boundary;
   vector<vector<double> > *dtrk_boundary;
   vector<vector<double> > *dvtx_boundary;
   vector<vector<double> > *dshr_boundary;
   vector<vector<double> > *dmc_boundary;
   Float_t         CosmicIP;
   Float_t         CosmicIPAll3D;
   Float_t         CosmicDirAll3D;
   Float_t         CosmicIPAll2DEnds;
   Float_t         CosmicDirAll2DEnds;
   Float_t         CosmicIPAll2DOvlp;
   Float_t         CosmicDirAll2DOvlp;
   Float_t         leeweight;
   Float_t         true_pt;
   Float_t         true_pt_visible;
   Float_t         true_p;
   Float_t         true_p_visible;
   Float_t         true_e_visible;
   Float_t         _opfilter_pe_beam;
   Float_t         _opfilter_pe_veto;
   Int_t           nu_pdg;
   Int_t           ccnc;
   Int_t           nu_parent_pdg;
   Int_t           nu_hadron_pdg;
   Int_t           nu_decay_mode;
   Double_t        par_decay_vx;
   Double_t        par_decay_vy;
   Double_t        par_decay_vz;
   Double_t        par_decay_px;
   Double_t        par_decay_py;
   Double_t        par_decay_pz;
   Double_t        baseline;
   Float_t         true_nu_px;
   Float_t         true_nu_py;
   Float_t         true_nu_pz;
   Int_t           swtrig_pre;
   Int_t           swtrig_post;
   Int_t           swtrig_pre_ext;
   Int_t           swtrig_post_ext;
   Int_t           interaction;
   Float_t         nu_e;
   Float_t         nu_l;
   Float_t         nu_pt;
   Float_t         theta;
   Bool_t          isVtxInFiducial;
   Bool_t          truthFiducial;
   Float_t         true_nu_vtx_t;
   Float_t         true_nu_vtx_x;
   Float_t         true_nu_vtx_y;
   Float_t         true_nu_vtx_z;
   Float_t         true_nu_vtx_sce_x;
   Float_t         true_nu_vtx_sce_y;
   Float_t         true_nu_vtx_sce_z;
   Float_t         reco_nu_vtx_x;
   Float_t         reco_nu_vtx_y;
   Float_t         reco_nu_vtx_z;
   Float_t         reco_nu_vtx_sce_x;
   Float_t         reco_nu_vtx_sce_y;
   Float_t         reco_nu_vtx_sce_z;
   Int_t           nmuon;
   Float_t         muon_e;
   Float_t         muon_c;
   Float_t         muon_p;
   Int_t           nelec;
   Float_t         elec_e;
   Float_t         elec_c;
   Float_t         elec_p;
   Float_t         elec_vx;
   Float_t         elec_vy;
   Float_t         elec_vz;
   Float_t         elec_px;
   Float_t         elec_py;
   Float_t         elec_pz;
   Int_t           npi0;
   Float_t         pi0_e;
   Float_t         pi0_c;
   Float_t         pi0_p;
   Int_t           nneutron;
   Int_t           nproton;
   Float_t         proton_e;
   Float_t         proton_c;
   Float_t         proton_p;
   Int_t           npion;
   Float_t         pion_e;
   Float_t         pion_c;
   Float_t         pion_p;
   Int_t           neta;
   Float_t         eta_e;
   Int_t           nslice;
   Int_t           crtveto;
   Float_t         crthitpe;
   vector<int>     *pfp_slice_idx;
   Int_t           category;
   vector<int>     *backtracked_pdg;
   vector<float>   *backtracked_e;
   vector<int>     *backtracked_tid;
   vector<float>   *backtracked_purity;
   vector<float>   *backtracked_completeness;
   vector<float>   *backtracked_overlay_purity;
   vector<float>   *backtracked_px;
   vector<float>   *backtracked_py;
   vector<float>   *backtracked_pz;
   vector<float>   *backtracked_start_x;
   vector<float>   *backtracked_start_y;
   vector<float>   *backtracked_start_z;
   vector<float>   *backtracked_start_t;
   vector<float>   *backtracked_start_U;
   vector<float>   *backtracked_start_V;
   vector<float>   *backtracked_start_Y;
   vector<float>   *backtracked_sce_start_x;
   vector<float>   *backtracked_sce_start_y;
   vector<float>   *backtracked_sce_start_z;
   vector<float>   *backtracked_sce_start_U;
   vector<float>   *backtracked_sce_start_V;
   vector<float>   *backtracked_sce_start_Y;
   Float_t         lep_e;
   Int_t           pass;
   Int_t           swtrig;
   Int_t           evnhits;
   Int_t           slpdg;
   Int_t           slnhits;
   Int_t           n_pfps;
   Int_t           n_tracks;
   Int_t           n_showers;
   vector<unsigned int> *pfp_generation_v;
   vector<unsigned int> *pfp_trk_daughters_v;
   vector<unsigned int> *pfp_shr_daughters_v;
   vector<unsigned int> *pfp_n_descendents_v;
   vector<float>   *pfp_vtx_x_v;
   vector<float>   *pfp_vtx_y_v;
   vector<float>   *pfp_vtx_z_v;
   vector<float>   *trk_score_v;
   vector<int>     *pfpdg;
   vector<int>     *pfnhits;
   vector<int>     *pfnplanehits_U;
   vector<int>     *pfnplanehits_V;
   vector<int>     *pfnplanehits_Y;
   vector<int>     *pfpplanesubclusters_U;
   vector<int>     *pfpplanesubclusters_V;
   vector<int>     *pfpplanesubclusters_Y;
   vector<float>   *pfpplanesubhitfracmax_U;
   vector<float>   *pfpplanesubhitfracmax_V;
   vector<float>   *pfpplanesubhitfracmax_Y;
   UInt_t          hits_u;
   UInt_t          hits_v;
   UInt_t          hits_y;
   UInt_t          slice_id;
   vector<float>   *slice_topo_score_v;
   Float_t         topological_score;
   Float_t         slclustfrac;
   vector<int>     *mc_pdg;
   vector<float>   *mc_E;
   vector<unsigned int> *mc_n_elastic;
   vector<unsigned int> *mc_n_inelastic;
   vector<float>   *mc_vx;
   vector<float>   *mc_vy;
   vector<float>   *mc_vz;
   vector<float>   *mc_endx;
   vector<float>   *mc_endy;
   vector<float>   *mc_endz;
   vector<float>   *mc_px;
   vector<float>   *mc_py;
   vector<float>   *mc_pz;
   vector<float>   *mc_end_p;
   vector<float>   *mc_completeness;
   vector<float>   *mc_purity;
   string          *endmuonprocess;
   Float_t         endmuonmichel;
   Int_t           filter_antibdt;
   Int_t           filter_ncpi0;
   Int_t           filter_pi0;
   Int_t           filter_ccinclusive;
   map<string,vector<double> > *weights;
   vector<unsigned short> *weightsFlux;
   vector<unsigned short> *weightsGenie;
   vector<unsigned short> *weightsReint;
   Float_t         weightSpline;
   Float_t         weightTune;
   Float_t         weightSplineTimesTune;
   Float_t         ppfx_cv;
   vector<unsigned short> *weightsPPFX;
   Double_t        knobRPAup;
   Double_t        knobRPAdn;
   Double_t        knobCCMECup;
   Double_t        knobCCMECdn;
   Double_t        knobAxFFCCQEup;
   Double_t        knobAxFFCCQEdn;
   Double_t        knobVecFFCCQEup;
   Double_t        knobVecFFCCQEdn;
   Double_t        knobDecayAngMECup;
   Double_t        knobDecayAngMECdn;
   Double_t        knobThetaDelta2Npiup;
   Double_t        knobThetaDelta2Npidn;
   Double_t        knobThetaDelta2NRadup;
   Double_t        knobThetaDelta2NRaddn;
   Double_t        knobNormCCCOHup;
   Double_t        knobNormCCCOHdn;
   Double_t        knobNormNCCOHup;
   Double_t        knobNormNCCOHdn;
   Double_t        knobxsr_scc_Fv3up;
   Double_t        knobxsr_scc_Fv3dn;
   Double_t        knobxsr_scc_Fa3up;
   Double_t        knobxsr_scc_Fa3dn;
   Double_t        RootinoFix;
   Float_t         flash_pe_flash_matching;
   vector<float>   *flash_pe_flash_matching_v;
   vector<float>   *slice_pe_flash_matching_v;
   Float_t         flash_time_flash_matching;
   Float_t         flash_y_flash_matching;
   Float_t         flash_z_flash_matching;
   Float_t         flash_timewidth_flash_matching;
   Float_t         flash_ywidth_flash_matching;
   Float_t         flash_zwidth_flash_matching;
   Float_t         nu_flashmatch_score;
   Float_t         nu_centerX;
   Float_t         nu_centerY;
   Float_t         nu_centerZ;
   Float_t         nu_totalCharge;
   Float_t         best_cosmic_flashmatch_score;
   Float_t         best_obviouscosmic_flashmatch_score;
   vector<float>   *cosmic_flashmatch_score_v;
   vector<float>   *cosmic_topological_score_v;
   vector<float>   *cosmic_centerX_v;
   vector<float>   *cosmic_centerY_v;
   vector<float>   *cosmic_centerZ_v;
   vector<float>   *cosmic_totalCharge_v;
   vector<int>     *cosmic_nhits_v;
   vector<int>     *cosmic_nunhits_v;
   vector<int>     *cosmic_isclear_v;
   Float_t         mcf_nu_e;
   Float_t         mcf_lep_e;
   Int_t           mcf_actvol;
   Int_t           mcf_nmm;
   Int_t           mcf_nmp;
   Int_t           mcf_nem;
   Int_t           mcf_nep;
   Int_t           mcf_np0;
   Int_t           mcf_npp;
   Int_t           mcf_npm;
   Float_t         mcf_mcshr_elec_etot;
   Int_t           mcf_pass_ccpi0;
   Int_t           mcf_pass_ncpi0;
   Int_t           mcf_pass_ccnopi;
   Int_t           mcf_pass_ncnopi;
   Int_t           mcf_pass_cccpi;
   Int_t           mcf_pass_nccpi;
   vector<float>   *X_SpcPts_v;
   vector<float>   *Y_SpcPts_v;
   vector<float>   *Z_SpcPts_v;
   UInt_t          shr_id_MCStool;
   UInt_t          shr_hits_max_MCStool;
   UInt_t          n_showers_contained_MCStool;
   vector<float>   *trkshrscore_v;
   Float_t         shrPCA_1Cr;
   Float_t         shrPCA_2Cr;
   Float_t         shrPCA_3Cr;
   Float_t         shrPCA_1Ce;
   Float_t         shrPCA_2Ce;
   Float_t         shrPCA_3Ce;
   Float_t         shrPCA1CAS;
   Float_t         shrPCA2CAS;
   Float_t         shrPCA3CAS;
   Float_t         shrPCA_1Cr2h;
   Float_t         shrPCA_2Cr2h;
   Float_t         shrPCA_3Cr2h;
   Float_t         shrPCA_1Cr1h;
   Float_t         shrPCA_2Cr1h;
   Float_t         shrPCA_3Cr1h;
   Float_t         shrMCSMom;
   Float_t         shrMCSMom1h;
   Float_t         shrMCSMom2h;
   Float_t         shrPCALen;
   UInt_t          n_shrSpcPts;
   vector<float>   *PCAWin_1Cr_5cm;
   vector<float>   *PCAWin_2Cr_5cm;
   vector<float>   *PCAWin_3Cr_5cm;
   vector<float>   *PCAWin_dist_5cm;
   vector<int>     *PCAWin_npts_5cm;
   Float_t         shrStart_5cm;
   Float_t         shrStartMCS_5cm;
   Float_t         shrMCSAS_5cm;
   Float_t         shrPCA1CAS_5cm;
   Float_t         shrPCA2CAS_5cm;
   Float_t         shrPCA3CAS_5cm;
   Float_t         shrPCA1CMed_5cm;
   vector<float>   *PCAWin_1Cr_2_5cm;
   vector<float>   *PCAWin_2Cr_2_5cm;
   vector<float>   *PCAWin_3Cr_2_5cm;
   vector<float>   *PCAWin_dist_2_5cm;
   vector<int>     *PCAWin_npts_2_5cm;
   Float_t         shrStart_2_5cm;
   Float_t         shrStartMCS_2_5cm;
   Float_t         shrMCSAS_2_5cm;
   Float_t         shrPCA1CAS_2_5cm;
   Float_t         shrPCA2CAS_2_5cm;
   Float_t         shrPCA3CAS_2_5cm;
   Float_t         shrPCA1CMed_2_5cm;
   Float_t         DeltaMed;
   Float_t         DeltaMed1h;
   Float_t         DeltaMed2h;
   Float_t         DeltaRMS;
   Float_t         DeltaRMS1h;
   Float_t         DeltaRMS2h;
   Float_t         CylFrac_1cm;
   Float_t         CylFrac1h_1cm;
   Float_t         CylFrac2h_1cm;
   Float_t         CylFrac_2cm;
   Float_t         CylFrac1h_2cm;
   Float_t         CylFrac2h_2cm;
   Float_t         CylFrac_3cm;
   Float_t         CylFrac1h_3cm;
   Float_t         CylFrac2h_3cm;
   Float_t         CylFrac_4cm;
   Float_t         CylFrac1h_4cm;
   Float_t         CylFrac2h_4cm;
   Float_t         CylFrac_5cm;
   Float_t         CylFrac1h_5cm;
   Float_t         CylFrac2h_5cm;
   Float_t         NeutrinoEnergy0;
   Float_t         NeutrinoEnergy1;
   Float_t         NeutrinoEnergy2;
   Float_t         SliceCaloEnergy0;
   Float_t         SliceCaloEnergy1;
   Float_t         SliceCaloEnergy2;
   Float_t         pi0_mcgamma0_e;
   Float_t         pi0_mcgamma0_px;
   Float_t         pi0_mcgamma0_py;
   Float_t         pi0_mcgamma0_pz;
   Float_t         pi0_mcrcdot0;
   Float_t         pi0_mcrce0;
   Float_t         pi0_mcgamma1_e;
   Float_t         pi0_mcgamma1_px;
   Float_t         pi0_mcgamma1_py;
   Float_t         pi0_mcgamma1_pz;
   Float_t         pi0_mcrcdot1;
   Float_t         pi0_mcrce1;
   Int_t           pi0_nshower;
   Int_t           pi0_ntrack;
   Int_t           pi0_ngamma;
   Float_t         pi0_radlen1;
   Float_t         pi0_radlen2;
   Float_t         pi0_dot1;
   Float_t         pi0_dot2;
   Float_t         pi0_energy1_Y;
   Float_t         pi0_energy2_Y;
   Float_t         pi0_dir1_x;
   Float_t         pi0_dir1_y;
   Float_t         pi0_dir1_z;
   Float_t         pi0_dir2_x;
   Float_t         pi0_dir2_y;
   Float_t         pi0_dir2_z;
   Float_t         pi0_dedx1_Y;
   Float_t         pi0_dedx2_Y;
   Float_t         pi0_dedx1_fit_Y;
   Float_t         pi0_dedx2_fit_Y;
   Float_t         pi0_energy1_V;
   Float_t         pi0_energy2_V;
   Float_t         pi0_dedx1_V;
   Float_t         pi0_dedx2_V;
   Float_t         pi0_dedx1_fit_V;
   Float_t         pi0_dedx2_fit_V;
   Float_t         pi0_energy1_U;
   Float_t         pi0_energy2_U;
   Float_t         pi0_dedx1_U;
   Float_t         pi0_dedx2_U;
   Float_t         pi0_dedx1_fit_U;
   Float_t         pi0_dedx2_fit_U;
   Float_t         pi0_shrscore1;
   Float_t         pi0_shrscore2;
   Float_t         pi0_gammadot;
   Float_t         pi0_mass_Y;
   Float_t         pi0_mass_V;
   Float_t         pi0_mass_U;
   Float_t         pi0_rc_vtx_x;
   Float_t         pi0_rc_vtx_y;
   Float_t         pi0_rc_vtx_z;
   Int_t           pi0truth_gamma_parent;
   Float_t         pi0truth_elec_edep;
   Float_t         pi0truth_elec_etot;
   Float_t         pi0truth_elec_dist;
   Int_t           pi0truth_elec_parent;
   Int_t           pi0truth_gamma1_tid;
   Float_t         pi0truth_gamma1_edep;
   Float_t         pi0truth_gamma1_etot;
   Float_t         pi0truth_gamma1_dist;
   Float_t         pi0truth_gamma1_elec1;
   Float_t         pi0truth_gamma1_elec2;
   Float_t         pi0truth_gamma1_xpos;
   Float_t         pi0truth_gamma1_ypos;
   Float_t         pi0truth_gamma1_zpos;
   Int_t           pi0truth_gamma2_tid;
   Float_t         pi0truth_gamma2_edep;
   Float_t         pi0truth_gamma2_etot;
   Float_t         pi0truth_gamma2_dist;
   Float_t         pi0truth_gamma2_elec1;
   Float_t         pi0truth_gamma2_elec2;
   Float_t         pi0truth_gamma2_xpos;
   Float_t         pi0truth_gamma2_ypos;
   Float_t         pi0truth_gamma2_zpos;
   Float_t         pi0truth_gammadot;
   Int_t           pi0truth_run;
   Int_t           pi0truth_sub;
   Int_t           pi0truth_evt;
   Int_t           nflag_pl1;
   Int_t           nnoise_pl1;
   Int_t           nslhits_pl1;
   Int_t           nslnoise_pl1;
   Int_t           nhits_pl1;
   Float_t         frac_slnoise_pl1;
   Float_t         secondshower_U_charge;
   Int_t           secondshower_U_nhit;
   Float_t         secondshower_U_vtxdist;
   Float_t         secondshower_U_eigenratio;
   Float_t         secondshower_U_dot;
   Float_t         secondshower_U_dir;
   Float_t         secondshower_V_charge;
   Int_t           secondshower_V_nhit;
   Float_t         secondshower_V_vtxdist;
   Float_t         secondshower_V_eigenratio;
   Float_t         secondshower_V_dot;
   Float_t         secondshower_V_dir;
   Float_t         secondshower_Y_charge;
   Int_t           secondshower_Y_nhit;
   Float_t         secondshower_Y_vtxdist;
   Float_t         secondshower_Y_eigenratio;
   Float_t         secondshower_Y_dot;
   Float_t         secondshower_Y_dir;
   vector<float>   *shr_dedx_u_v;
   vector<float>   *shr_dedx_v_v;
   vector<float>   *shr_dedx_y_v;
   vector<float>   *shr_energy_u_v;
   vector<float>   *shr_energy_v_v;
   vector<float>   *shr_energy_y_v;
   vector<unsigned long> *shr_pfp_id_v;
   vector<float>   *shr_start_x_v;
   vector<float>   *shr_start_y_v;
   vector<float>   *shr_start_z_v;
   vector<float>   *shr_dist_v;
   vector<float>   *shr_start_U_v;
   vector<float>   *shr_start_V_v;
   vector<float>   *shr_px_v;
   vector<float>   *shr_py_v;
   vector<float>   *shr_pz_v;
   vector<float>   *shr_openangle_v;
   vector<float>   *shr_theta_v;
   vector<float>   *shr_phi_v;
   vector<float>   *shr_pitch_u_v;
   vector<float>   *shr_pitch_v_v;
   vector<float>   *shr_pitch_y_v;
   vector<int>     *shr_tkfit_nhits_v;
   vector<float>   *shr_tkfit_start_x_v;
   vector<float>   *shr_tkfit_start_y_v;
   vector<float>   *shr_tkfit_start_z_v;
   vector<float>   *shr_tkfit_start_U_v;
   vector<float>   *shr_tkfit_start_V_v;
   vector<float>   *shr_tkfit_theta_v;
   vector<float>   *shr_tkfit_phi_v;
   vector<float>   *shr_tkfit_pitch_u_v;
   vector<float>   *shr_tkfit_pitch_v_v;
   vector<float>   *shr_tkfit_pitch_y_v;
   vector<float>   *shr_tkfit_dedx_u_v;
   vector<float>   *shr_tkfit_dedx_v_v;
   vector<float>   *shr_tkfit_dedx_y_v;
   vector<float>   *shr_tkfit_gap10_dedx_u_v;
   vector<float>   *shr_tkfit_gap10_dedx_v_v;
   vector<float>   *shr_tkfit_gap10_dedx_y_v;
   vector<int>     *shr_tkfit_dedx_nhits_u_v;
   vector<int>     *shr_tkfit_dedx_nhits_v_v;
   vector<int>     *shr_tkfit_dedx_nhits_y_v;
   vector<float>   *shr_llr_pid_u_v;
   vector<float>   *shr_llr_pid_v_v;
   vector<float>   *shr_llr_pid_y_v;
   vector<float>   *shr_llr_pid_v;
   vector<float>   *shr_llr_pid_score_v;
   vector<float>   *shr_moliere_avg_v;
   vector<float>   *shr_moliere_rms_v;
   Int_t           origevnunhits;
   Int_t           evnunhits;
   Int_t           evlepnhits;
   Int_t           evpronhits;
   Int_t           evpi1nhits;
   Int_t           evpi0nhits;
   Int_t           evneunhits;
   Int_t           evgamnhits;
   Int_t           evothnhits;
   Int_t           slnunhits;
   Int_t           sllepnhits;
   Int_t           slpronhits;
   Int_t           slpi1nhits;
   Int_t           slpi0nhits;
   Int_t           slneunhits;
   Int_t           slgamnhits;
   Int_t           slothnhits;
   vector<int>     *pfnunhits;
   vector<int>     *pflepnhits;
   vector<int>     *pfpronhits;
   vector<int>     *pfpi1nhits;
   vector<int>     *pfpi0nhits;
   vector<int>     *pfneunhits;
   vector<int>     *pfgamnhits;
   vector<int>     *pfothnhits;
   Float_t         nu_completeness_from_pfp;
   Float_t         nu_purity_from_pfp;
   vector<float>   *trk_bragg_p_v;
   vector<float>   *trk_bragg_mu_v;
   vector<float>   *trk_bragg_pion_v;
   vector<float>   *trk_bragg_mip_v;
   vector<float>   *trk_bragg_p_alt_dir_v;
   vector<float>   *trk_bragg_mu_alt_dir_v;
   vector<float>   *trk_bragg_pion_alt_dir_v;
   vector<bool>    *trk_bragg_p_fwd_preferred_v;
   vector<bool>    *trk_bragg_mu_fwd_preferred_v;
   vector<bool>    *trk_bragg_pion_fwd_preferred_v;
   vector<float>   *trk_pida_v;
   vector<float>   *trk_pid_chipr_v;
   vector<float>   *trk_pid_chipi_v;
   vector<float>   *trk_pid_chika_v;
   vector<float>   *trk_pid_chimu_v;
   vector<float>   *trk_bragg_p_u_v;
   vector<float>   *trk_bragg_mu_u_v;
   vector<float>   *trk_bragg_pion_u_v;
   vector<float>   *trk_bragg_mip_u_v;
   vector<float>   *trk_bragg_p_alt_dir_u_v;
   vector<float>   *trk_bragg_mu_alt_dir_u_v;
   vector<float>   *trk_bragg_pion_alt_dir_u_v;
   vector<bool>    *trk_bragg_p_fwd_preferred_u_v;
   vector<bool>    *trk_bragg_mu_fwd_preferred_u_v;
   vector<bool>    *trk_bragg_pion_fwd_preferred_u_v;
   vector<float>   *trk_pida_u_v;
   vector<float>   *trk_pid_chipr_u_v;
   vector<float>   *trk_pid_chipi_u_v;
   vector<float>   *trk_pid_chika_u_v;
   vector<float>   *trk_pid_chimu_u_v;
   vector<float>   *trk_bragg_p_v_v;
   vector<float>   *trk_bragg_mu_v_v;
   vector<float>   *trk_bragg_pion_v_v;
   vector<float>   *trk_bragg_mip_v_v;
   vector<float>   *trk_bragg_p_alt_dir_v_v;
   vector<float>   *trk_bragg_mu_alt_dir_v_v;
   vector<float>   *trk_bragg_pion_alt_dir_v_v;
   vector<bool>    *trk_bragg_p_fwd_preferred_v_v;
   vector<bool>    *trk_bragg_mu_fwd_preferred_v_v;
   vector<bool>    *trk_bragg_pion_fwd_preferred_v_v;
   vector<float>   *trk_pida_v_v;
   vector<float>   *trk_pid_chipr_v_v;
   vector<float>   *trk_pid_chipi_v_v;
   vector<float>   *trk_pid_chika_v_v;
   vector<float>   *trk_pid_chimu_v_v;
   vector<unsigned long> *trk_pfp_id_v;
   vector<float>   *trk_dir_x_v;
   vector<float>   *trk_dir_y_v;
   vector<float>   *trk_dir_z_v;
   vector<float>   *trk_start_x_v;
   vector<float>   *trk_start_y_v;
   vector<float>   *trk_start_z_v;
   vector<float>   *trk_sce_start_x_v;
   vector<float>   *trk_sce_start_y_v;
   vector<float>   *trk_sce_start_z_v;
   vector<float>   *trk_end_x_v;
   vector<float>   *trk_end_y_v;
   vector<float>   *trk_end_z_v;
   vector<float>   *trk_sce_end_x_v;
   vector<float>   *trk_sce_end_y_v;
   vector<float>   *trk_sce_end_z_v;
   vector<float>   *trk_distance_v;
   vector<float>   *trk_theta_v;
   vector<float>   *trk_phi_v;
   vector<float>   *trk_len_v;
   vector<float>   *trk_mcs_muon_mom_v;
   vector<float>   *trk_range_muon_mom_v;
   vector<float>   *trk_energy_proton_v;
   vector<float>   *trk_energy_muon_v;
   vector<float>   *trk_calo_energy_u_v;
   vector<float>   *trk_calo_energy_v_v;
   vector<float>   *trk_calo_energy_y_v;
   vector<float>   *trk_llr_pid_u_v;
   vector<float>   *trk_llr_pid_v_v;
   vector<float>   *trk_llr_pid_y_v;
   vector<float>   *trk_llr_pid_v;
   vector<float>   *trk_llr_pid_score_v;
   vector<float>   *trk_trunk_dEdx_u_v;
   vector<float>   *trk_trunk_dEdx_v_v;
   vector<float>   *trk_trunk_dEdx_y_v;
   vector<float>   *trk_trunk_rr_dEdx_u_v;
   vector<float>   *trk_trunk_rr_dEdx_v_v;
   vector<float>   *trk_trunk_rr_dEdx_y_v;
   vector<int>     *trk_nhits_u_v;
   vector<int>     *trk_nhits_v_v;
   vector<int>     *trk_nhits_y_v;
   vector<float>   *trk_avg_deflection_mean_v;
   vector<float>   *trk_avg_deflection_stdev_v;
   vector<float>   *trk_avg_deflection_separation_mean_v;
   vector<int>     *trk_end_spacepoints_v;
   Float_t         bdt_nuNCpi0;
   Float_t         bdt_numuCCpi0;
   Float_t         bdt_numuCC;
   Float_t         bdt_ext;
   Float_t         bdt_cosmic;
   Float_t         bdt_global;
   Int_t           pass_antibdt_filter;
   Float_t         bdt_pi0_np;
   Float_t         bdt_nonpi0_np;
   Float_t         bdt_bkg_0p;
   Float_t         anglediff_Y;
   Float_t         anglediff_V;
   Float_t         anglediff_U;
   Float_t         trkpid;

   // List of branches
   TBranch        *b_selected;   //!
   TBranch        *b_run;   //!
   TBranch        *b_sub;   //!
   TBranch        *b_evt;   //!
   TBranch        *b_trk_pfp_id;   //!
   TBranch        *b_shr_pfp_id;   //!
   TBranch        *b_trk2_pfp_id;   //!
   TBranch        *b_shr2_pfp_id;   //!
   TBranch        *b_trk3_pfp_id;   //!
   TBranch        *b_shr3_pfp_id;   //!
   TBranch        *b_shr_energy_tot;   //!
   TBranch        *b_shr_energy;   //!
   TBranch        *b_shr_energy_second;   //!
   TBranch        *b_shr_energy_third;   //!
   TBranch        *b_shr_energy_tot_cali;   //!
   TBranch        *b_shr_energy_cali;   //!
   TBranch        *b_shr_energy_second_cali;   //!
   TBranch        *b_shr_energy_third_cali;   //!
   TBranch        *b_shr_theta;   //!
   TBranch        *b_shr_phi;   //!
   TBranch        *b_shr_pca_0;   //!
   TBranch        *b_shr_pca_1;   //!
   TBranch        *b_shr_pca_2;   //!
   TBranch        *b_shr_px;   //!
   TBranch        *b_shr_py;   //!
   TBranch        *b_shr_pz;   //!
   TBranch        *b_shr_openangle;   //!
   TBranch        *b_shr_tkfit_start_x;   //!
   TBranch        *b_shr_tkfit_start_y;   //!
   TBranch        *b_shr_tkfit_start_z;   //!
   TBranch        *b_shr_tkfit_theta;   //!
   TBranch        *b_shr_tkfit_phi;   //!
   TBranch        *b_shr_start_x;   //!
   TBranch        *b_shr_start_y;   //!
   TBranch        *b_shr_start_z;   //!
   TBranch        *b_shr_dedx_Y;   //!
   TBranch        *b_shr_dedx_V;   //!
   TBranch        *b_shr_dedx_U;   //!
   TBranch        *b_shr_dedx_Y_cali;   //!
   TBranch        *b_shr_dedx_V_cali;   //!
   TBranch        *b_shr_dedx_U_cali;   //!
   TBranch        *b_shr_tkfit_dedx_Y;   //!
   TBranch        *b_shr_tkfit_dedx_V;   //!
   TBranch        *b_shr_tkfit_dedx_U;   //!
   TBranch        *b_shr_tkfit_dedx_max;   //!
   TBranch        *b_shr_tkfit_nhits_Y;   //!
   TBranch        *b_shr_tkfit_nhits_V;   //!
   TBranch        *b_shr_tkfit_nhits_U;   //!
   TBranch        *b_shr_llrpid_dedx_Y;   //!
   TBranch        *b_shr_llrpid_dedx_V;   //!
   TBranch        *b_shr_llrpid_dedx_U;   //!
   TBranch        *b_shr_llrpid_dedx;   //!
   TBranch        *b_shr_tkfit_dedx_Y_alt;   //!
   TBranch        *b_shr_tkfit_dedx_V_alt;   //!
   TBranch        *b_shr_tkfit_dedx_U_alt;   //!
   TBranch        *b_shr_tkfit_nhits_Y_alt;   //!
   TBranch        *b_shr_tkfit_nhits_V_alt;   //!
   TBranch        *b_shr_tkfit_nhits_U_alt;   //!
   TBranch        *b__trkfit;   //!
   TBranch        *b_shr_tkfit_npoints;   //!
   TBranch        *b_shr_tkfit_npointsvalid;   //!
   TBranch        *b_shr_trkfitmedangle;   //!
   TBranch        *b_shrmoliereavg;   //!
   TBranch        *b_shrmoliererms;   //!
   TBranch        *b_shr1shr2moliereavg;   //!
   TBranch        *b_shr1shr2moliererms;   //!
   TBranch        *b_shr1trk1moliereavg;   //!
   TBranch        *b_shr1trk1moliererms;   //!
   TBranch        *b_shr1trk2moliereavg;   //!
   TBranch        *b_shr1trk2moliererms;   //!
   TBranch        *b_ismerged;   //!
   TBranch        *b_merge_bestdot;   //!
   TBranch        *b_merge_bestdist;   //!
   TBranch        *b_merge_vtx_x;   //!
   TBranch        *b_merge_vtx_y;   //!
   TBranch        *b_merge_vtx_z;   //!
   TBranch        *b_merge_tk_ipfp;   //!
   TBranch        *b_shr_tkfit_2cm_dedx_Y;   //!
   TBranch        *b_shr_tkfit_2cm_dedx_V;   //!
   TBranch        *b_shr_tkfit_2cm_dedx_U;   //!
   TBranch        *b_shr_tkfit_2cm_nhits_Y;   //!
   TBranch        *b_shr_tkfit_2cm_nhits_V;   //!
   TBranch        *b_shr_tkfit_2cm_nhits_U;   //!
   TBranch        *b_shr_tkfit_gap05_dedx_Y;   //!
   TBranch        *b_shr_tkfit_gap05_dedx_V;   //!
   TBranch        *b_shr_tkfit_gap05_dedx_U;   //!
   TBranch        *b_shr_tkfit_gap05_nhits_Y;   //!
   TBranch        *b_shr_tkfit_gap05_nhits_V;   //!
   TBranch        *b_shr_tkfit_gap05_nhits_U;   //!
   TBranch        *b_shr_tkfit_gap10_dedx_Y;   //!
   TBranch        *b_shr_tkfit_gap10_dedx_V;   //!
   TBranch        *b_shr_tkfit_gap10_dedx_U;   //!
   TBranch        *b_shr_tkfit_gap10_nhits_Y;   //!
   TBranch        *b_shr_tkfit_gap10_nhits_V;   //!
   TBranch        *b_shr_tkfit_gap10_nhits_U;   //!
   TBranch        *b_shr_chipr;   //!
   TBranch        *b_shr_chimu;   //!
   TBranch        *b_shr_bragg_p;   //!
   TBranch        *b_shr_bragg_mu;   //!
   TBranch        *b_shr_bragg_mip;   //!
   TBranch        *b_shr_bragg_kaon;   //!
   TBranch        *b_shr_bragg_pion;   //!
   TBranch        *b_tksh_distance;   //!
   TBranch        *b_tksh_angle;   //!
   TBranch        *b_shr_distance;   //!
   TBranch        *b_shr_score;   //!
   TBranch        *b_shr_bkt_pdg;   //!
   TBranch        *b_shr_bkt_purity;   //!
   TBranch        *b_shr_bkt_completeness;   //!
   TBranch        *b_shr_bkt_E;   //!
   TBranch        *b_trk_len;   //!
   TBranch        *b_trk_theta;   //!
   TBranch        *b_trk_phi;   //!
   TBranch        *b_trk_energy;   //!
   TBranch        *b_trk_energy_muon;   //!
   TBranch        *b_trk_energy_muon_mcs;   //!
   TBranch        *b_trk_energy_tot;   //!
   TBranch        *b_trk_energy_muon_tot;   //!
   TBranch        *b_trk_distance;   //!
   TBranch        *b_trk_score;   //!
   TBranch        *b_trk_bkt_pdg;   //!
   TBranch        *b_trk_bkt_purity;   //!
   TBranch        *b_trk_bkt_completeness;   //!
   TBranch        *b_trk_bkt_E;   //!
   TBranch        *b_trk_chipr_best;   //!
   TBranch        *b_trk_chipr_worst;   //!
   TBranch        *b_trk_chimu_best;   //!
   TBranch        *b_trk_chimu_worst;   //!
   TBranch        *b_trk_chipr;   //!
   TBranch        *b_trk_chimu;   //!
   TBranch        *b_trk_pida;   //!
   TBranch        *b_trk_bragg_p;   //!
   TBranch        *b_trk_bragg_mu;   //!
   TBranch        *b_trk_bragg_mip;   //!
   TBranch        *b_trk_bragg_kaon;   //!
   TBranch        *b_trk_bragg_pion;   //!
   TBranch        *b_trk_hits_max;   //!
   TBranch        *b_shr_hits_max;   //!
   TBranch        *b_all_shr_hits;   //!
   TBranch        *b_all_trk_hits;   //!
   TBranch        *b_all_shr_energies;   //!
   TBranch        *b_all_trk_energies;   //!
   TBranch        *b_shr_hits_max;   //!
   TBranch        *b_trk_hits_2nd;   //!
   TBranch        *b_shr_hits_2nd;   //!
   TBranch        *b_trk_hits_3rd;   //!
   TBranch        *b_shr_hits_3rd;   //!
   TBranch        *b_trkshrhitdist0;   //!
   TBranch        *b_trkshrhitdist1;   //!
   TBranch        *b_trkshrhitdist2;   //!
   TBranch        *b_trk2shrhitdist0;   //!
   TBranch        *b_trk2shrhitdist1;   //!
   TBranch        *b_trk2shrhitdist2;   //!
   TBranch        *b_trk1trk2hitdist0;   //!
   TBranch        *b_trk1trk2hitdist1;   //!
   TBranch        *b_trk1trk2hitdist2;   //!
   TBranch        *b_total_hits_y;   //!
   TBranch        *b_extra_energy_y;   //!
   TBranch        *b_trk_energy_hits_tot;   //!
   TBranch        *b_subcluster;   //!
   TBranch        *b_shrsubclusters0;   //!
   TBranch        *b_shrsubclusters1;   //!
   TBranch        *b_shrsubclusters2;   //!
   TBranch        *b_shrclusfrac0;   //!
   TBranch        *b_shrclusfrac1;   //!
   TBranch        *b_shrclusfrac2;   //!
   TBranch        *b_shrclusdir0;   //!
   TBranch        *b_shrclusdir1;   //!
   TBranch        *b_shrclusdir2;   //!
   TBranch        *b_shr_hits_tot;   //!
   TBranch        *b_shr_hits_y_tot;   //!
   TBranch        *b_shr_hits_u_tot;   //!
   TBranch        *b_shr_hits_v_tot;   //!
   TBranch        *b_trk_hits_tot;   //!
   TBranch        *b_trk_hits_y_tot;   //!
   TBranch        *b_trk_hits_u_tot;   //!
   TBranch        *b_trk_hits_v_tot;   //!
   TBranch        *b_elecclusters_U_charge;   //!
   TBranch        *b_elecclusters_V_charge;   //!
   TBranch        *b_elecclusters_Y_charge;   //!
   TBranch        *b_elecclusters_U_N;   //!
   TBranch        *b_elecclusters_V_N;   //!
   TBranch        *b_elecclusters_Y_N;   //!
   TBranch        *b_n_tracks_contained;   //!
   TBranch        *b_n_showers_contained;   //!
   TBranch        *b_matched_E;   //!
   TBranch        *b_hits_ratio;   //!
   TBranch        *b_contained_fraction;   //!
   TBranch        *b_sps_contained_fraction;   //!
   TBranch        *b_pt;   //!
   TBranch        *b_p;   //!
   TBranch        *b_pt_assume_muon;   //!
   TBranch        *b_p_assume_muon;   //!
   TBranch        *b_reco_e;   //!
   TBranch        *b_dvtx;   //!
   TBranch        *b_dtrk;   //!
   TBranch        *b_contained_sps_ratio;   //!
   TBranch        *b_dtrk_x_boundary;   //!
   TBranch        *b_dtrk_y_boundary;   //!
   TBranch        *b_dtrk_z_boundary;   //!
   TBranch        *b_dshr_x_boundary;   //!
   TBranch        *b_dshr_y_boundary;   //!
   TBranch        *b_dshr_z_boundary;   //!
   TBranch        *b_dvtx_x_boundary;   //!
   TBranch        *b_dvtx_y_boundary;   //!
   TBranch        *b_dvtx_z_boundary;   //!
   TBranch        *b_dtrk_boundary;   //!
   TBranch        *b_dvtx_boundary;   //!
   TBranch        *b_dshr_boundary;   //!
   TBranch        *b_dmc_boundary;   //!
   TBranch        *b_CosmicIP;   //!
   TBranch        *b_CosmicIPAll3D;   //!
   TBranch        *b_CosmicDirAll3D;   //!
   TBranch        *b_CosmicIPAll2DEnds;   //!
   TBranch        *b_CosmicDirAll2DEnds;   //!
   TBranch        *b_CosmicIPAll2DOvlp;   //!
   TBranch        *b_CosmicDirAll2DOvlp;   //!
   TBranch        *b_leeweight;   //!
   TBranch        *b_true_pt;   //!
   TBranch        *b_true_pt_visible;   //!
   TBranch        *b_true_p;   //!
   TBranch        *b_true_p_visible;   //!
   TBranch        *b_true_e_visible;   //!
   TBranch        *b_opfilter_pe_beam;   //!
   TBranch        *b_opfilter_pe_veto;   //!
   TBranch        *b_nu_pdg;   //!
   TBranch        *b_ccnc;   //!
   TBranch        *b_nu_parent_pdg;   //!
   TBranch        *b_nu_hadron_pdg;   //!
   TBranch        *b_nu_decay_mode;   //!
   TBranch        *b_par_decay_vx;   //!
   TBranch        *b_par_decay_vy;   //!
   TBranch        *b_par_decay_vz;   //!
   TBranch        *b_par_decay_px;   //!
   TBranch        *b_par_decay_py;   //!
   TBranch        *b_par_decay_pz;   //!
   TBranch        *b_baseline;   //!
   TBranch        *b_true_nu_px;   //!
   TBranch        *b_true_nu_py;   //!
   TBranch        *b_true_nu_pz;   //!
   TBranch        *b_swtrig_pre;   //!
   TBranch        *b_swtrig_post;   //!
   TBranch        *b_swtrig_pre_ext;   //!
   TBranch        *b_swtrig_post_ext;   //!
   TBranch        *b_interaction;   //!
   TBranch        *b_nu_e;   //!
   TBranch        *b_nu_l;   //!
   TBranch        *b_nu_pt;   //!
   TBranch        *b_theta;   //!
   TBranch        *b_isVtxInFiducial;   //!
   TBranch        *b_truthFiducial;   //!
   TBranch        *b_true_nu_vtx_t;   //!
   TBranch        *b_true_nu_vtx_x;   //!
   TBranch        *b_true_nu_vtx_y;   //!
   TBranch        *b_true_nu_vtx_z;   //!
   TBranch        *b_true_nu_vtx_sce_x;   //!
   TBranch        *b_true_nu_vtx_sce_y;   //!
   TBranch        *b_true_nu_vtx_sce_z;   //!
   TBranch        *b_reco_nu_vtx_x;   //!
   TBranch        *b_reco_nu_vtx_y;   //!
   TBranch        *b_reco_nu_vtx_z;   //!
   TBranch        *b_reco_nu_vtx_sce_x;   //!
   TBranch        *b_reco_nu_vtx_sce_y;   //!
   TBranch        *b_reco_nu_vtx_sce_z;   //!
   TBranch        *b_nmuon;   //!
   TBranch        *b_muon_e;   //!
   TBranch        *b_muon_c;   //!
   TBranch        *b_muon_p;   //!
   TBranch        *b_nelec;   //!
   TBranch        *b_elec_e;   //!
   TBranch        *b_elec_c;   //!
   TBranch        *b_elec_p;   //!
   TBranch        *b_elec_vx;   //!
   TBranch        *b_elec_vy;   //!
   TBranch        *b_elec_vz;   //!
   TBranch        *b_elec_px;   //!
   TBranch        *b_elec_py;   //!
   TBranch        *b_elec_pz;   //!
   TBranch        *b_npi0;   //!
   TBranch        *b_pi0_e;   //!
   TBranch        *b_pi0_c;   //!
   TBranch        *b_pi0_p;   //!
   TBranch        *b_nneutron;   //!
   TBranch        *b_nproton;   //!
   TBranch        *b_proton_e;   //!
   TBranch        *b_proton_c;   //!
   TBranch        *b_proton_p;   //!
   TBranch        *b_npion;   //!
   TBranch        *b_pion_e;   //!
   TBranch        *b_pion_c;   //!
   TBranch        *b_pion_p;   //!
   TBranch        *b_neta;   //!
   TBranch        *b_eta_e;   //!
   TBranch        *b_nslice;   //!
   TBranch        *b_crtveto;   //!
   TBranch        *b_crthitpe;   //!
   TBranch        *b_pfp_slice_idx;   //!
   TBranch        *b_category;   //!
   TBranch        *b_backtracked_pdg;   //!
   TBranch        *b_backtracked_e;   //!
   TBranch        *b_backtracked_tid;   //!
   TBranch        *b_backtracked_purity;   //!
   TBranch        *b_backtracked_completeness;   //!
   TBranch        *b_backtracked_overlay_purity;   //!
   TBranch        *b_backtracked_px;   //!
   TBranch        *b_backtracked_py;   //!
   TBranch        *b_backtracked_pz;   //!
   TBranch        *b_backtracked_start_x;   //!
   TBranch        *b_backtracked_start_y;   //!
   TBranch        *b_backtracked_start_z;   //!
   TBranch        *b_backtracked_start_t;   //!
   TBranch        *b_backtracked_start_U;   //!
   TBranch        *b_backtracked_start_V;   //!
   TBranch        *b_backtracked_start_Y;   //!
   TBranch        *b_backtracked_sce_start_x;   //!
   TBranch        *b_backtracked_sce_start_y;   //!
   TBranch        *b_backtracked_sce_start_z;   //!
   TBranch        *b_backtracked_sce_start_U;   //!
   TBranch        *b_backtracked_sce_start_V;   //!
   TBranch        *b_backtracked_sce_start_Y;   //!
   TBranch        *b_lep_e;   //!
   TBranch        *b_pass;   //!
   TBranch        *b_swtrig;   //!
   TBranch        *b_evnhits;   //!
   TBranch        *b_slpdg;   //!
   TBranch        *b_slnhits;   //!
   TBranch        *b_n_pfps;   //!
   TBranch        *b_n_tracks;   //!
   TBranch        *b_n_showers;   //!
   TBranch        *b_pfp_generation_v;   //!
   TBranch        *b_pfp_trk_daughters_v;   //!
   TBranch        *b_pfp_shr_daughters_v;   //!
   TBranch        *b_pfp_n_descendents_v;   //!
   TBranch        *b_pfp_vtx_x_v;   //!
   TBranch        *b_pfp_vtx_y_v;   //!
   TBranch        *b_pfp_vtx_z_v;   //!
   TBranch        *b_trk_score_v;   //!
   TBranch        *b_pfpdg;   //!
   TBranch        *b_pfnhits;   //!
   TBranch        *b_pfnplanehits_U;   //!
   TBranch        *b_pfnplanehits_V;   //!
   TBranch        *b_pfnplanehits_Y;   //!
   TBranch        *b_pfpplanesubclusters_U;   //!
   TBranch        *b_pfpplanesubclusters_V;   //!
   TBranch        *b_pfpplanesubclusters_Y;   //!
   TBranch        *b_pfpplanesubhitfracmax_U;   //!
   TBranch        *b_pfpplanesubhitfracmax_V;   //!
   TBranch        *b_pfpplanesubhitfracmax_Y;   //!
   TBranch        *b_hits_u;   //!
   TBranch        *b_hits_v;   //!
   TBranch        *b_hits_y;   //!
   TBranch        *b_slice_id;   //!
   TBranch        *b_slice_topo_score_v;   //!
   TBranch        *b_topological_score;   //!
   TBranch        *b_slclustfrac;   //!
   TBranch        *b_mc_pdg;   //!
   TBranch        *b_mc_E;   //!
   TBranch        *b_mc_n_elastic;   //!
   TBranch        *b_mc_n_inelastic;   //!
   TBranch        *b_mc_vx;   //!
   TBranch        *b_mc_vy;   //!
   TBranch        *b_mc_vz;   //!
   TBranch        *b_mc_endx;   //!
   TBranch        *b_mc_endy;   //!
   TBranch        *b_mc_endz;   //!
   TBranch        *b_mc_px;   //!
   TBranch        *b_mc_py;   //!
   TBranch        *b_mc_pz;   //!
   TBranch        *b_mc_end_p;   //!
   TBranch        *b_mc_completeness;   //!
   TBranch        *b_mc_purity;   //!
   TBranch        *b_endmuonprocess;   //!
   TBranch        *b_endmuonmichel;   //!
   TBranch        *b_filter_antibdt;   //!
   TBranch        *b_filter_ncpi0;   //!
   TBranch        *b_filter_pi0;   //!
   TBranch        *b_filter_ccinclusive;   //!
   TBranch        *b_weights;   //!
   TBranch        *b_weightsFlux;   //!
   TBranch        *b_weightsGenie;   //!
   TBranch        *b_weightsReint;   //!
   TBranch        *b_weightSpline;   //!
   TBranch        *b_weightTune;   //!
   TBranch        *b_weightSplineTimesTune;   //!
   TBranch        *b_ppfx_cv;   //!
   TBranch        *b_weightsPPFX;   //!
   TBranch        *b_knobRPAup;   //!
   TBranch        *b_knobRPAdn;   //!
   TBranch        *b_knobCCMECup;   //!
   TBranch        *b_knobCCMECdn;   //!
   TBranch        *b_knobAxFFCCQEup;   //!
   TBranch        *b_knobAxFFCCQEdn;   //!
   TBranch        *b_knobVecFFCCQEup;   //!
   TBranch        *b_knobVecFFCCQEdn;   //!
   TBranch        *b_knobDecayAngMECup;   //!
   TBranch        *b_knobDecayAngMECdn;   //!
   TBranch        *b_knobThetaDelta2Npiup;   //!
   TBranch        *b_knobThetaDelta2Npidn;   //!
   TBranch        *b_knobThetaDelta2NRadup;   //!
   TBranch        *b_knobThetaDelta2NRaddn;   //!
   TBranch        *b_knobNormCCCOHup;   //!
   TBranch        *b_knobNormCCCOHdn;   //!
   TBranch        *b_knobNormNCCOHup;   //!
   TBranch        *b_knobNormNCCOHdn;   //!
   TBranch        *b_knobxsr_scc_Fv3up;   //!
   TBranch        *b_knobxsr_scc_Fv3dn;   //!
   TBranch        *b_knobxsr_scc_Fa3up;   //!
   TBranch        *b_knobxsr_scc_Fa3dn;   //!
   TBranch        *b_RootinoFix;   //!
   TBranch        *b_flash_pe_flash_matching;   //!
   TBranch        *b_flash_pe_flash_matching_v;   //!
   TBranch        *b_slice_pe_flash_matching_v;   //!
   TBranch        *b_flash_time_flash_matching;   //!
   TBranch        *b_flash_y_flash_matching;   //!
   TBranch        *b_flash_z_flash_matching;   //!
   TBranch        *b_flash_timewidth_flash_matching;   //!
   TBranch        *b_flash_ywidth_flash_matching;   //!
   TBranch        *b_flash_zwidth_flash_matching;   //!
   TBranch        *b_nu_flashmatch_score;   //!
   TBranch        *b_nu_centerX;   //!
   TBranch        *b_nu_centerY;   //!
   TBranch        *b_nu_centerZ;   //!
   TBranch        *b_nu_totalCharge;   //!
   TBranch        *b_best_cosmic_flashmatch_score;   //!
   TBranch        *b_best_obviouscosmic_flashmatch_score;   //!
   TBranch        *b_cosmic_flashmatch_score_v;   //!
   TBranch        *b_cosmic_topological_score_v;   //!
   TBranch        *b_cosmic_centerX_v;   //!
   TBranch        *b_cosmic_centerY_v;   //!
   TBranch        *b_cosmic_centerZ_v;   //!
   TBranch        *b_cosmic_totalCharge_v;   //!
   TBranch        *b_cosmic_nhits_v;   //!
   TBranch        *b_cosmic_nunhits_v;   //!
   TBranch        *b_cosmic_isclear_v;   //!
   TBranch        *b_mcf_nu_e;   //!
   TBranch        *b_mcf_lep_e;   //!
   TBranch        *b_mcf_actvol;   //!
   TBranch        *b_mcf_nmm;   //!
   TBranch        *b_mcf_nmp;   //!
   TBranch        *b_mcf_nem;   //!
   TBranch        *b_mcf_nep;   //!
   TBranch        *b_mcf_np0;   //!
   TBranch        *b_mcf_npp;   //!
   TBranch        *b_mcf_npm;   //!
   TBranch        *b_mcf_mcshr_elec_etot;   //!
   TBranch        *b_mcf_pass_ccpi0;   //!
   TBranch        *b_mcf_pass_ncpi0;   //!
   TBranch        *b_mcf_pass_ccnopi;   //!
   TBranch        *b_mcf_pass_ncnopi;   //!
   TBranch        *b_mcf_pass_cccpi;   //!
   TBranch        *b_mcf_pass_nccpi;   //!
   TBranch        *b_X_SpcPts_v;   //!
   TBranch        *b_Y_SpcPts_v;   //!
   TBranch        *b_Z_SpcPts_v;   //!
   TBranch        *b_shr_pfp_id;   //!
   TBranch        *b_shr_hits_max_MCStool;   //!
   TBranch        *b_n_showers_contained_MCStool;   //!
   TBranch        *b_trkshrscore_v;   //!
   TBranch        *b_shrPCA_1Cr;   //!
   TBranch        *b_shrPCA_2Cr;   //!
   TBranch        *b_shrPCA_3Cr;   //!
   TBranch        *b_shrPCA_1Ce;   //!
   TBranch        *b_shrPCA_2Ce;   //!
   TBranch        *b_shrPCA_3Ce;   //!
   TBranch        *b_shrPCA1CAS;   //!
   TBranch        *b_shrPCA2CAS;   //!
   TBranch        *b_shrPCA3CAS;   //!
   TBranch        *b_shrPCA_1Cr2h;   //!
   TBranch        *b_shrPCA_2Cr2h;   //!
   TBranch        *b_shrPCA_3Cr2h;   //!
   TBranch        *b_shrPCA_1Cr1h;   //!
   TBranch        *b_shrPCA_2Cr1h;   //!
   TBranch        *b_shrPCA_3Cr1h;   //!
   TBranch        *b_shrMCSMom;   //!
   TBranch        *b_shrMCSMom1h;   //!
   TBranch        *b_shrMCSMom2h;   //!
   TBranch        *b_shrPCALen;   //!
   TBranch        *b_n_shrSpcPts;   //!
   TBranch        *b_PCAWin_1Cr_5cm;   //!
   TBranch        *b_PCAWin_2Cr_5cm;   //!
   TBranch        *b_PCAWin_3Cr_5cm;   //!
   TBranch        *b_PCAWin_dist_5cm;   //!
   TBranch        *b_PCAWin_npts_5cm;   //!
   TBranch        *b_shrStart_5cm;   //!
   TBranch        *b_shrStartMCS_5cm;   //!
   TBranch        *b_shrMCSAS_5cm;   //!
   TBranch        *b_shrPCA1CAS_5cm;   //!
   TBranch        *b_shrPCA2CAS_5cm;   //!
   TBranch        *b_shrPCA3CAS_5cm;   //!
   TBranch        *b__shrPCA1CMed_5cm;   //!
   TBranch        *b_PCAWin_1Cr_2_5cm;   //!
   TBranch        *b_PCAWin_2Cr_2_5cm;   //!
   TBranch        *b_PCAWin_3Cr_2_5cm;   //!
   TBranch        *b_PCAWin_dist_2_5cm;   //!
   TBranch        *b_PCAWin_npts_2_5cm;   //!
   TBranch        *b_shrStart_2_5cm;   //!
   TBranch        *b_shrStartMCS_2_5cm;   //!
   TBranch        *b_shrMCSAS_2_5cm;   //!
   TBranch        *b_shrPCA1CAS_2_5cm;   //!
   TBranch        *b_shrPCA2CAS_2_5cm;   //!
   TBranch        *b_shrPCA3CAS_2_5cm;   //!
   TBranch        *b__shrPCA1CMed_2_5cm;   //!
   TBranch        *b_DeltaMed;   //!
   TBranch        *b_DeltaMed1h;   //!
   TBranch        *b_DeltaMed2h;   //!
   TBranch        *b_DeltaRMS;   //!
   TBranch        *b_DeltaRMS1h;   //!
   TBranch        *b_DeltaRMS2h;   //!
   TBranch        *b_CylFrac_1cm;   //!
   TBranch        *b_CylFrac1h_1cm;   //!
   TBranch        *b_CylFrac2h_1cm;   //!
   TBranch        *b_CylFrac_2cm;   //!
   TBranch        *b_CylFrac1h_2cm;   //!
   TBranch        *b_CylFrac2h_2cm;   //!
   TBranch        *b_CylFrac_3cm;   //!
   TBranch        *b_CylFrac1h_3cm;   //!
   TBranch        *b_CylFrac2h_3cm;   //!
   TBranch        *b_CylFrac_4cm;   //!
   TBranch        *b_CylFrac1h_4cm;   //!
   TBranch        *b_CylFrac2h_4cm;   //!
   TBranch        *b_CylFrac_5cm;   //!
   TBranch        *b_CylFrac1h_5cm;   //!
   TBranch        *b_CylFrac2h_5cm;   //!
   TBranch        *b_NeutrinoEnergy0;   //!
   TBranch        *b_NeutrinoEnergy1;   //!
   TBranch        *b_NeutrinoEnergy2;   //!
   TBranch        *b_SliceCaloEnergy0;   //!
   TBranch        *b_SliceCaloEnergy1;   //!
   TBranch        *b_SliceCaloEnergy2;   //!
   TBranch        *b_pi0_mcgamma0_e;   //!
   TBranch        *b_pi0_mcgamma0_px;   //!
   TBranch        *b_pi0_mcgamma0_py;   //!
   TBranch        *b_pi0_mcgamma0_pz;   //!
   TBranch        *b_pi0_mcrcdot0;   //!
   TBranch        *b_pi0_mcrce0;   //!
   TBranch        *b_pi0_mcgamma1_e;   //!
   TBranch        *b_pi0_mcgamma1_px;   //!
   TBranch        *b_pi0_mcgamma1_py;   //!
   TBranch        *b_pi0_mcgamma1_pz;   //!
   TBranch        *b_pi0_mcrcdot1;   //!
   TBranch        *b_pi0_mcrce1;   //!
   TBranch        *b_pi0_nshower;   //!
   TBranch        *b_pi0_ntrack;   //!
   TBranch        *b_pi0_ngamma;   //!
   TBranch        *b_pi0_radlen1;   //!
   TBranch        *b_pi0_radlen2;   //!
   TBranch        *b_pi0_dot1;   //!
   TBranch        *b_pi0_dot2;   //!
   TBranch        *b_pi0_energy1_Y;   //!
   TBranch        *b_pi0_energy2_Y;   //!
   TBranch        *b_pi0_dir1_x;   //!
   TBranch        *b_pi0_dir1_y;   //!
   TBranch        *b_pi0_dir1_z;   //!
   TBranch        *b_pi0_dir2_x;   //!
   TBranch        *b_pi0_dir2_y;   //!
   TBranch        *b_pi0_dir2_z;   //!
   TBranch        *b_pi0_dedx1_Y;   //!
   TBranch        *b_pi0_dedx2_Y;   //!
   TBranch        *b_pi0_dedx1_fit_Y;   //!
   TBranch        *b_pi0_dedx2_fit_Y;   //!
   TBranch        *b_pi0_energy1_V;   //!
   TBranch        *b_pi0_energy2_V;   //!
   TBranch        *b_pi0_dedx1_V;   //!
   TBranch        *b_pi0_dedx2_V;   //!
   TBranch        *b_pi0_dedx1_fit_V;   //!
   TBranch        *b_pi0_dedx2_fit_V;   //!
   TBranch        *b_pi0_energy1_U;   //!
   TBranch        *b_pi0_energy2_U;   //!
   TBranch        *b_pi0_dedx1_U;   //!
   TBranch        *b_pi0_dedx2_U;   //!
   TBranch        *b_pi0_dedx1_fit_U;   //!
   TBranch        *b_pi0_dedx2_fit_U;   //!
   TBranch        *b_pi0_shrscore1;   //!
   TBranch        *b_pi0_shrscore2;   //!
   TBranch        *b_pi0_gammadot;   //!
   TBranch        *b_pi0_mass_Y;   //!
   TBranch        *b_pi0_mass_V;   //!
   TBranch        *b_pi0_mass_U;   //!
   TBranch        *b_pi0_rc_vtx_x;   //!
   TBranch        *b_pi0_rc_vtx_y;   //!
   TBranch        *b_pi0_rc_vtx_z;   //!
   TBranch        *b_pi0truth_gamma_parent;   //!
   TBranch        *b_pi0truth_elec_edep;   //!
   TBranch        *b_pi0truth_elec_etot;   //!
   TBranch        *b_pi0truth_elec_dist;   //!
   TBranch        *b_pi0truth_elec_parent;   //!
   TBranch        *b_pi0truth_gamma1_tid;   //!
   TBranch        *b_pi0truth_gamma1_edep;   //!
   TBranch        *b_pi0truth_gamma1_etot;   //!
   TBranch        *b_pi0truth_gamma1_dist;   //!
   TBranch        *b_pi0truth_gamma1_elec1;   //!
   TBranch        *b_pi0truth_gamma1_elec2;   //!
   TBranch        *b_pi0truth_gamma1_xpos;   //!
   TBranch        *b_pi0truth_gamma1_ypos;   //!
   TBranch        *b_pi0truth_gamma1_zpos;   //!
   TBranch        *b_pi0truth_gamma2_tid;   //!
   TBranch        *b_pi0truth_gamma2_edep;   //!
   TBranch        *b_pi0truth_gamma2_etot;   //!
   TBranch        *b_pi0truth_gamma2_dist;   //!
   TBranch        *b_pi0truth_gamma2_elec1;   //!
   TBranch        *b_pi0truth_gamma2_elec2;   //!
   TBranch        *b_pi0truth_gamma2_xpos;   //!
   TBranch        *b_pi0truth_gamma2_ypos;   //!
   TBranch        *b_pi0truth_gamma2_zpos;   //!
   TBranch        *b_pi0truth_gammadot;   //!
   TBranch        *b_pi0truth_run;   //!
   TBranch        *b_pi0truth_sub;   //!
   TBranch        *b_pi0truth_evt;   //!
   TBranch        *b_nflag_pl1;   //!
   TBranch        *b_nnoise_pl1;   //!
   TBranch        *b_nslhits_pl1;   //!
   TBranch        *b_nslnoise_pl1;   //!
   TBranch        *b_nhits_pl1;   //!
   TBranch        *b_frac_slnoise_pl1;   //!
   TBranch        *b_secondshower_U_charge;   //!
   TBranch        *b_secondshower_U_nhit;   //!
   TBranch        *b_secondshower_U_vtxdist;   //!
   TBranch        *b_secondshower_U_eigenratio;   //!
   TBranch        *b_secondshower_U_dot;   //!
   TBranch        *b_secondshower_U_dir;   //!
   TBranch        *b_secondshower_V_charge;   //!
   TBranch        *b_secondshower_V_nhit;   //!
   TBranch        *b_secondshower_V_vtxdist;   //!
   TBranch        *b_secondshower_V_eigenratio;   //!
   TBranch        *b_secondshower_V_dot;   //!
   TBranch        *b_secondshower_V_dir;   //!
   TBranch        *b_secondshower_Y_charge;   //!
   TBranch        *b_secondshower_Y_nhit;   //!
   TBranch        *b_secondshower_Y_vtxdist;   //!
   TBranch        *b_secondshower_Y_eigenratio;   //!
   TBranch        *b_secondshower_Y_dot;   //!
   TBranch        *b_secondshower_Y_dir;   //!
   TBranch        *b_shr_dedx_u_v;   //!
   TBranch        *b_shr_dedx_v_v;   //!
   TBranch        *b_shr_dedx_y_v;   //!
   TBranch        *b_shr_energy_u_v;   //!
   TBranch        *b_shr_energy_v_v;   //!
   TBranch        *b_shr_energy_y_v;   //!
   TBranch        *b_shr_pfp_id_v;   //!
   TBranch        *b_shr_start_x_v;   //!
   TBranch        *b_shr_start_y_v;   //!
   TBranch        *b_shr_start_z_v;   //!
   TBranch        *b_shr_dist_v;   //!
   TBranch        *b_shr_start_U_v;   //!
   TBranch        *b_shr_start_V_v;   //!
   TBranch        *b_shr_px_v;   //!
   TBranch        *b_shr_py_v;   //!
   TBranch        *b_shr_pz_v;   //!
   TBranch        *b_shr_openangle_v;   //!
   TBranch        *b_shr_theta_v;   //!
   TBranch        *b_shr_phi_v;   //!
   TBranch        *b_shr_pitch_u_v;   //!
   TBranch        *b_shr_pitch_v_v;   //!
   TBranch        *b_shr_pitch_y_v;   //!
   TBranch        *b_shr_tkfit_nhits_v;   //!
   TBranch        *b_shr_tkfit_start_x_v;   //!
   TBranch        *b_shr_tkfit_start_y_v;   //!
   TBranch        *b_shr_tkfit_start_z_v;   //!
   TBranch        *b_shr_tkfit_start_U_v;   //!
   TBranch        *b_shr_tkfit_start_V_v;   //!
   TBranch        *b_shr_tkfit_theta_v;   //!
   TBranch        *b_shr_tkfit_phi_v;   //!
   TBranch        *b_shr_tkfit_pitch_u_v;   //!
   TBranch        *b_shr_tkfit_pitch_v_v;   //!
   TBranch        *b_shr_tkfit_pitch_y_v;   //!
   TBranch        *b_shr_tkfit_dedx_u_v;   //!
   TBranch        *b_shr_tkfit_dedx_v_v;   //!
   TBranch        *b_shr_tkfit_dedx_y_v;   //!
   TBranch        *b_shr_tkfit_gap10_dedx_u_v;   //!
   TBranch        *b_shr_tkfit_gap10_dedx_v_v;   //!
   TBranch        *b_shr_tkfit_gap10_dedx_y_v;   //!
   TBranch        *b_shr_tkfit_dedx_nhits_u_v;   //!
   TBranch        *b_shr_tkfit_dedx_nhits_v_v;   //!
   TBranch        *b_shr_tkfit_dedx_nhits_y_v;   //!
   TBranch        *b_shr_llr_pid_u_v;   //!
   TBranch        *b_shr_llr_pid_v_v;   //!
   TBranch        *b_shr_llr_pid_y_v;   //!
   TBranch        *b_shr_llr_pid_v;   //!
   TBranch        *b_shr_llr_pid_score_v;   //!
   TBranch        *b_shr_moliere_avg_v;   //!
   TBranch        *b_shr_moliere_rms_v;   //!
   TBranch        *b_origevnunhits;   //!
   TBranch        *b_evnunhits;   //!
   TBranch        *b_evlepnhits;   //!
   TBranch        *b_evpronhits;   //!
   TBranch        *b_evpi1nhits;   //!
   TBranch        *b_evpi0nhits;   //!
   TBranch        *b_evneunhits;   //!
   TBranch        *b_evgamnhits;   //!
   TBranch        *b_evothnhits;   //!
   TBranch        *b_slnunhits;   //!
   TBranch        *b_sllepnhits;   //!
   TBranch        *b_slpronhits;   //!
   TBranch        *b_slpi1nhits;   //!
   TBranch        *b_slpi0nhits;   //!
   TBranch        *b_slneunhits;   //!
   TBranch        *b_slgamnhits;   //!
   TBranch        *b_slothnhits;   //!
   TBranch        *b_pfnunhits;   //!
   TBranch        *b_pflepnhits;   //!
   TBranch        *b_pfpronhits;   //!
   TBranch        *b_pfpi1nhits;   //!
   TBranch        *b_pfpi0nhits;   //!
   TBranch        *b_pfneunhits;   //!
   TBranch        *b_pfgamnhits;   //!
   TBranch        *b_pfothnhits;   //!
   TBranch        *b_nu_completeness_from_pfp;   //!
   TBranch        *b_nu_purity_from_pfp;   //!
   TBranch        *b_trk_bragg_p_v;   //!
   TBranch        *b_trk_bragg_mu_v;   //!
   TBranch        *b_trk_bragg_pion_v;   //!
   TBranch        *b_trk_bragg_mip_v;   //!
   TBranch        *b_trk_bragg_p_alt_dir_v;   //!
   TBranch        *b_trk_bragg_mu_alt_dir_v;   //!
   TBranch        *b_trk_bragg_pion_alt_dir_v;   //!
   TBranch        *b_trk_bragg_p_fwd_preferred_v;   //!
   TBranch        *b_trk_bragg_mu_fwd_preferred_v;   //!
   TBranch        *b_trk_bragg_pion_fwd_preferred_v;   //!
   TBranch        *b_trk_pida_v;   //!
   TBranch        *b_trk_pid_chipr_v;   //!
   TBranch        *b_trk_pid_chipi_v;   //!
   TBranch        *b_trk_pid_chika_v;   //!
   TBranch        *b_trk_pid_chimu_v;   //!
   TBranch        *b_trk_bragg_p_u_v;   //!
   TBranch        *b_trk_bragg_mu_u_v;   //!
   TBranch        *b_trk_bragg_pion_u_v;   //!
   TBranch        *b_trk_bragg_mip_u_v;   //!
   TBranch        *b_trk_bragg_p_alt_dir_u_v;   //!
   TBranch        *b_trk_bragg_mu_alt_dir_u_v;   //!
   TBranch        *b_trk_bragg_pion_alt_dir_u_v;   //!
   TBranch        *b_trk_bragg_p_fwd_preferred_u_v;   //!
   TBranch        *b_trk_bragg_mu_fwd_preferred_u_v;   //!
   TBranch        *b_trk_bragg_pion_fwd_preferred_u_v;   //!
   TBranch        *b_trk_pida_u_v;   //!
   TBranch        *b_trk_pid_chipr_u_v;   //!
   TBranch        *b_trk_pid_chipi_u_v;   //!
   TBranch        *b_trk_pid_chika_u_v;   //!
   TBranch        *b_trk_pid_chimu_u_v;   //!
   TBranch        *b_trk_bragg_p_v_v;   //!
   TBranch        *b_trk_bragg_mu_v_v;   //!
   TBranch        *b_trk_bragg_pion_v_v;   //!
   TBranch        *b_trk_bragg_mip_v_v;   //!
   TBranch        *b_trk_bragg_p_alt_dir_v_v;   //!
   TBranch        *b_trk_bragg_mu_alt_dir_v_v;   //!
   TBranch        *b_trk_bragg_pion_alt_dir_v_v;   //!
   TBranch        *b_trk_bragg_p_fwd_preferred_v_v;   //!
   TBranch        *b_trk_bragg_mu_fwd_preferred_v_v;   //!
   TBranch        *b_trk_bragg_pion_fwd_preferred_v_v;   //!
   TBranch        *b_trk_pida_v_v;   //!
   TBranch        *b_trk_pid_chipr_v_v;   //!
   TBranch        *b_trk_pid_chipi_v_v;   //!
   TBranch        *b_trk_pid_chika_v_v;   //!
   TBranch        *b_trk_pid_chimu_v_v;   //!
   TBranch        *b_trk_pfp_id_v;   //!
   TBranch        *b_trk_dir_x_v;   //!
   TBranch        *b_trk_dir_y_v;   //!
   TBranch        *b_trk_dir_z_v;   //!
   TBranch        *b_trk_start_x_v;   //!
   TBranch        *b_trk_start_y_v;   //!
   TBranch        *b_trk_start_z_v;   //!
   TBranch        *b_trk_sce_start_x_v;   //!
   TBranch        *b_trk_sce_start_y_v;   //!
   TBranch        *b_trk_sce_start_z_v;   //!
   TBranch        *b_trk_end_x_v;   //!
   TBranch        *b_trk_end_y_v;   //!
   TBranch        *b_trk_end_z_v;   //!
   TBranch        *b_trk_sce_end_x_v;   //!
   TBranch        *b_trk_sce_end_y_v;   //!
   TBranch        *b_trk_sce_end_z_v;   //!
   TBranch        *b_trk_distance_v;   //!
   TBranch        *b_trk_theta_v;   //!
   TBranch        *b_trk_phi_v;   //!
   TBranch        *b_trk_len_v;   //!
   TBranch        *b_trk_mcs_muon_mom_v;   //!
   TBranch        *b_trk_range_muon_mom_v;   //!
   TBranch        *b_trk_energy_proton_v;   //!
   TBranch        *b_trk_energy_muon_v;   //!
   TBranch        *b_trk_calo_energy_u_v;   //!
   TBranch        *b_trk_calo_energy_v_v;   //!
   TBranch        *b_trk_calo_energy_y_v;   //!
   TBranch        *b_trk_llr_pid_u_v;   //!
   TBranch        *b_trk_llr_pid_v_v;   //!
   TBranch        *b_trk_llr_pid_y_v;   //!
   TBranch        *b_trk_llr_pid_v;   //!
   TBranch        *b_trk_llr_pid_score_v;   //!
   TBranch        *b_trk_trunk_dEdx_u_v;   //!
   TBranch        *b_trk_trunk_dEdx_v_v;   //!
   TBranch        *b_trk_trunk_dEdx_y_v;   //!
   TBranch        *b_trk_trunk_rr_dEdx_u_v;   //!
   TBranch        *b_trk_trunk_rr_dEdx_v_v;   //!
   TBranch        *b_trk_trunk_rr_dEdx_y_v;   //!
   TBranch        *b_trk_nhits_u_v;   //!
   TBranch        *b_trk_nhits_v_v;   //!
   TBranch        *b_trk_nhits_y_v;   //!
   TBranch        *b_trk_avg_deflection_mean_v;   //!
   TBranch        *b_trk_avg_deflection_stdev_v;   //!
   TBranch        *b_trk_avg_deflection_separation_mean_v;   //!
   TBranch        *b_trk_end_spacepoints_v;   //!
   TBranch        *b_bdt_nuNCpi0;   //!
   TBranch        *b_bdt_numuCCpi0;   //!
   TBranch        *b_bdt_numuCC;   //!
   TBranch        *b_bdt_ext;   //!
   TBranch        *b_bdt_cosmic;   //!
   TBranch        *b_bdt_global;   //!
   TBranch        *b_bdt_global;   //!
   TBranch        *b_bdt_pi0_np;   //!
   TBranch        *b_bdt_nonpi0_np;   //!
   TBranch        *b_bdt_bkg_0p;   //!
   TBranch        *b_anglediff_Y;   //!
   TBranch        *b_anglediff_V;   //!
   TBranch        *b_anglediff_U;   //!
   TBranch        *b_trkpid;   //!

   NeutrinoSelectionFilter(TTree *tree=0);
   virtual ~NeutrinoSelectionFilter();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual bool     Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef NeutrinoSelectionFilter_cxx
NeutrinoSelectionFilter::NeutrinoSelectionFilter(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Run1_fhc_new_numi_flux_fhc_pandora_ntuple.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("Run1_fhc_new_numi_flux_fhc_pandora_ntuple.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("Run1_fhc_new_numi_flux_fhc_pandora_ntuple.root:/nuselection");
      dir->GetObject("NeutrinoSelectionFilter",tree);

   }
   Init(tree);
}

NeutrinoSelectionFilter::~NeutrinoSelectionFilter()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t NeutrinoSelectionFilter::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t NeutrinoSelectionFilter::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void NeutrinoSelectionFilter::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   all_shr_hits = 0;
   all_trk_hits = 0;
   all_shr_energies = 0;
   all_trk_energies = 0;
   dtrk_x_boundary = 0;
   dtrk_y_boundary = 0;
   dtrk_z_boundary = 0;
   dshr_x_boundary = 0;
   dshr_y_boundary = 0;
   dshr_z_boundary = 0;
   dvtx_x_boundary = 0;
   dvtx_y_boundary = 0;
   dvtx_z_boundary = 0;
   dtrk_boundary = 0;
   dvtx_boundary = 0;
   dshr_boundary = 0;
   dmc_boundary = 0;
   pfp_slice_idx = 0;
   backtracked_pdg = 0;
   backtracked_e = 0;
   backtracked_tid = 0;
   backtracked_purity = 0;
   backtracked_completeness = 0;
   backtracked_overlay_purity = 0;
   backtracked_px = 0;
   backtracked_py = 0;
   backtracked_pz = 0;
   backtracked_start_x = 0;
   backtracked_start_y = 0;
   backtracked_start_z = 0;
   backtracked_start_t = 0;
   backtracked_start_U = 0;
   backtracked_start_V = 0;
   backtracked_start_Y = 0;
   backtracked_sce_start_x = 0;
   backtracked_sce_start_y = 0;
   backtracked_sce_start_z = 0;
   backtracked_sce_start_U = 0;
   backtracked_sce_start_V = 0;
   backtracked_sce_start_Y = 0;
   pfp_generation_v = 0;
   pfp_trk_daughters_v = 0;
   pfp_shr_daughters_v = 0;
   pfp_n_descendents_v = 0;
   pfp_vtx_x_v = 0;
   pfp_vtx_y_v = 0;
   pfp_vtx_z_v = 0;
   trk_score_v = 0;
   pfpdg = 0;
   pfnhits = 0;
   pfnplanehits_U = 0;
   pfnplanehits_V = 0;
   pfnplanehits_Y = 0;
   pfpplanesubclusters_U = 0;
   pfpplanesubclusters_V = 0;
   pfpplanesubclusters_Y = 0;
   pfpplanesubhitfracmax_U = 0;
   pfpplanesubhitfracmax_V = 0;
   pfpplanesubhitfracmax_Y = 0;
   slice_topo_score_v = 0;
   mc_pdg = 0;
   mc_E = 0;
   mc_n_elastic = 0;
   mc_n_inelastic = 0;
   mc_vx = 0;
   mc_vy = 0;
   mc_vz = 0;
   mc_endx = 0;
   mc_endy = 0;
   mc_endz = 0;
   mc_px = 0;
   mc_py = 0;
   mc_pz = 0;
   mc_end_p = 0;
   mc_completeness = 0;
   mc_purity = 0;
   endmuonprocess = 0;
   weights = 0;
   weightsFlux = 0;
   weightsGenie = 0;
   weightsReint = 0;
   weightsPPFX = 0;
   flash_pe_flash_matching_v = 0;
   slice_pe_flash_matching_v = 0;
   cosmic_flashmatch_score_v = 0;
   cosmic_topological_score_v = 0;
   cosmic_centerX_v = 0;
   cosmic_centerY_v = 0;
   cosmic_centerZ_v = 0;
   cosmic_totalCharge_v = 0;
   cosmic_nhits_v = 0;
   cosmic_nunhits_v = 0;
   cosmic_isclear_v = 0;
   X_SpcPts_v = 0;
   Y_SpcPts_v = 0;
   Z_SpcPts_v = 0;
   trkshrscore_v = 0;
   PCAWin_1Cr_5cm = 0;
   PCAWin_2Cr_5cm = 0;
   PCAWin_3Cr_5cm = 0;
   PCAWin_dist_5cm = 0;
   PCAWin_npts_5cm = 0;
   PCAWin_1Cr_2_5cm = 0;
   PCAWin_2Cr_2_5cm = 0;
   PCAWin_3Cr_2_5cm = 0;
   PCAWin_dist_2_5cm = 0;
   PCAWin_npts_2_5cm = 0;
   shr_dedx_u_v = 0;
   shr_dedx_v_v = 0;
   shr_dedx_y_v = 0;
   shr_energy_u_v = 0;
   shr_energy_v_v = 0;
   shr_energy_y_v = 0;
   shr_pfp_id_v = 0;
   shr_start_x_v = 0;
   shr_start_y_v = 0;
   shr_start_z_v = 0;
   shr_dist_v = 0;
   shr_start_U_v = 0;
   shr_start_V_v = 0;
   shr_px_v = 0;
   shr_py_v = 0;
   shr_pz_v = 0;
   shr_openangle_v = 0;
   shr_theta_v = 0;
   shr_phi_v = 0;
   shr_pitch_u_v = 0;
   shr_pitch_v_v = 0;
   shr_pitch_y_v = 0;
   shr_tkfit_nhits_v = 0;
   shr_tkfit_start_x_v = 0;
   shr_tkfit_start_y_v = 0;
   shr_tkfit_start_z_v = 0;
   shr_tkfit_start_U_v = 0;
   shr_tkfit_start_V_v = 0;
   shr_tkfit_theta_v = 0;
   shr_tkfit_phi_v = 0;
   shr_tkfit_pitch_u_v = 0;
   shr_tkfit_pitch_v_v = 0;
   shr_tkfit_pitch_y_v = 0;
   shr_tkfit_dedx_u_v = 0;
   shr_tkfit_dedx_v_v = 0;
   shr_tkfit_dedx_y_v = 0;
   shr_tkfit_gap10_dedx_u_v = 0;
   shr_tkfit_gap10_dedx_v_v = 0;
   shr_tkfit_gap10_dedx_y_v = 0;
   shr_tkfit_dedx_nhits_u_v = 0;
   shr_tkfit_dedx_nhits_v_v = 0;
   shr_tkfit_dedx_nhits_y_v = 0;
   shr_llr_pid_u_v = 0;
   shr_llr_pid_v_v = 0;
   shr_llr_pid_y_v = 0;
   shr_llr_pid_v = 0;
   shr_llr_pid_score_v = 0;
   shr_moliere_avg_v = 0;
   shr_moliere_rms_v = 0;
   pfnunhits = 0;
   pflepnhits = 0;
   pfpronhits = 0;
   pfpi1nhits = 0;
   pfpi0nhits = 0;
   pfneunhits = 0;
   pfgamnhits = 0;
   pfothnhits = 0;
   trk_bragg_p_v = 0;
   trk_bragg_mu_v = 0;
   trk_bragg_pion_v = 0;
   trk_bragg_mip_v = 0;
   trk_bragg_p_alt_dir_v = 0;
   trk_bragg_mu_alt_dir_v = 0;
   trk_bragg_pion_alt_dir_v = 0;
   trk_bragg_p_fwd_preferred_v = 0;
   trk_bragg_mu_fwd_preferred_v = 0;
   trk_bragg_pion_fwd_preferred_v = 0;
   trk_pida_v = 0;
   trk_pid_chipr_v = 0;
   trk_pid_chipi_v = 0;
   trk_pid_chika_v = 0;
   trk_pid_chimu_v = 0;
   trk_bragg_p_u_v = 0;
   trk_bragg_mu_u_v = 0;
   trk_bragg_pion_u_v = 0;
   trk_bragg_mip_u_v = 0;
   trk_bragg_p_alt_dir_u_v = 0;
   trk_bragg_mu_alt_dir_u_v = 0;
   trk_bragg_pion_alt_dir_u_v = 0;
   trk_bragg_p_fwd_preferred_u_v = 0;
   trk_bragg_mu_fwd_preferred_u_v = 0;
   trk_bragg_pion_fwd_preferred_u_v = 0;
   trk_pida_u_v = 0;
   trk_pid_chipr_u_v = 0;
   trk_pid_chipi_u_v = 0;
   trk_pid_chika_u_v = 0;
   trk_pid_chimu_u_v = 0;
   trk_bragg_p_v_v = 0;
   trk_bragg_mu_v_v = 0;
   trk_bragg_pion_v_v = 0;
   trk_bragg_mip_v_v = 0;
   trk_bragg_p_alt_dir_v_v = 0;
   trk_bragg_mu_alt_dir_v_v = 0;
   trk_bragg_pion_alt_dir_v_v = 0;
   trk_bragg_p_fwd_preferred_v_v = 0;
   trk_bragg_mu_fwd_preferred_v_v = 0;
   trk_bragg_pion_fwd_preferred_v_v = 0;
   trk_pida_v_v = 0;
   trk_pid_chipr_v_v = 0;
   trk_pid_chipi_v_v = 0;
   trk_pid_chika_v_v = 0;
   trk_pid_chimu_v_v = 0;
   trk_pfp_id_v = 0;
   trk_dir_x_v = 0;
   trk_dir_y_v = 0;
   trk_dir_z_v = 0;
   trk_start_x_v = 0;
   trk_start_y_v = 0;
   trk_start_z_v = 0;
   trk_sce_start_x_v = 0;
   trk_sce_start_y_v = 0;
   trk_sce_start_z_v = 0;
   trk_end_x_v = 0;
   trk_end_y_v = 0;
   trk_end_z_v = 0;
   trk_sce_end_x_v = 0;
   trk_sce_end_y_v = 0;
   trk_sce_end_z_v = 0;
   trk_distance_v = 0;
   trk_theta_v = 0;
   trk_phi_v = 0;
   trk_len_v = 0;
   trk_mcs_muon_mom_v = 0;
   trk_range_muon_mom_v = 0;
   trk_energy_proton_v = 0;
   trk_energy_muon_v = 0;
   trk_calo_energy_u_v = 0;
   trk_calo_energy_v_v = 0;
   trk_calo_energy_y_v = 0;
   trk_llr_pid_u_v = 0;
   trk_llr_pid_v_v = 0;
   trk_llr_pid_y_v = 0;
   trk_llr_pid_v = 0;
   trk_llr_pid_score_v = 0;
   trk_trunk_dEdx_u_v = 0;
   trk_trunk_dEdx_v_v = 0;
   trk_trunk_dEdx_y_v = 0;
   trk_trunk_rr_dEdx_u_v = 0;
   trk_trunk_rr_dEdx_v_v = 0;
   trk_trunk_rr_dEdx_y_v = 0;
   trk_nhits_u_v = 0;
   trk_nhits_v_v = 0;
   trk_nhits_y_v = 0;
   trk_avg_deflection_mean_v = 0;
   trk_avg_deflection_stdev_v = 0;
   trk_avg_deflection_separation_mean_v = 0;
   trk_end_spacepoints_v = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("selected", &selected, &b_selected);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("sub", &sub, &b_sub);
   fChain->SetBranchAddress("evt", &evt, &b_evt);
   fChain->SetBranchAddress("trk_id", &trk_id, &b_trk_pfp_id);
   fChain->SetBranchAddress("shr_id", &shr_id, &b_shr_pfp_id);
   fChain->SetBranchAddress("trk2_id", &trk2_id, &b_trk2_pfp_id);
   fChain->SetBranchAddress("shr2_id", &shr2_id, &b_shr2_pfp_id);
   fChain->SetBranchAddress("trk3_id", &trk3_id, &b_trk3_pfp_id);
   fChain->SetBranchAddress("shr3_id", &shr3_id, &b_shr3_pfp_id);
   fChain->SetBranchAddress("shr_energy_tot", &shr_energy_tot, &b_shr_energy_tot);
   fChain->SetBranchAddress("shr_energy", &shr_energy, &b_shr_energy);
   fChain->SetBranchAddress("shr_energy_second", &shr_energy_second, &b_shr_energy_second);
   fChain->SetBranchAddress("shr_energy_third", &shr_energy_third, &b_shr_energy_third);
   fChain->SetBranchAddress("shr_energy_tot_cali", &shr_energy_tot_cali, &b_shr_energy_tot_cali);
   fChain->SetBranchAddress("shr_energy_cali", &shr_energy_cali, &b_shr_energy_cali);
   fChain->SetBranchAddress("shr_energy_second_cali", &shr_energy_second_cali, &b_shr_energy_second_cali);
   fChain->SetBranchAddress("shr_energy_third_cali", &shr_energy_third_cali, &b_shr_energy_third_cali);
   fChain->SetBranchAddress("shr_theta", &shr_theta, &b_shr_theta);
   fChain->SetBranchAddress("shr_phi", &shr_phi, &b_shr_phi);
   fChain->SetBranchAddress("shr_pca_0", &shr_pca_0, &b_shr_pca_0);
   fChain->SetBranchAddress("shr_pca_1", &shr_pca_1, &b_shr_pca_1);
   fChain->SetBranchAddress("shr_pca_2", &shr_pca_2, &b_shr_pca_2);
   fChain->SetBranchAddress("shr_px", &shr_px, &b_shr_px);
   fChain->SetBranchAddress("shr_py", &shr_py, &b_shr_py);
   fChain->SetBranchAddress("shr_pz", &shr_pz, &b_shr_pz);
   fChain->SetBranchAddress("shr_openangle", &shr_openangle, &b_shr_openangle);
   fChain->SetBranchAddress("shr_tkfit_start_x", &shr_tkfit_start_x, &b_shr_tkfit_start_x);
   fChain->SetBranchAddress("shr_tkfit_start_y", &shr_tkfit_start_y, &b_shr_tkfit_start_y);
   fChain->SetBranchAddress("shr_tkfit_start_z", &shr_tkfit_start_z, &b_shr_tkfit_start_z);
   fChain->SetBranchAddress("shr_tkfit_theta", &shr_tkfit_theta, &b_shr_tkfit_theta);
   fChain->SetBranchAddress("shr_tkfit_phi", &shr_tkfit_phi, &b_shr_tkfit_phi);
   fChain->SetBranchAddress("shr_start_x", &shr_start_x, &b_shr_start_x);
   fChain->SetBranchAddress("shr_start_y", &shr_start_y, &b_shr_start_y);
   fChain->SetBranchAddress("shr_start_z", &shr_start_z, &b_shr_start_z);
   fChain->SetBranchAddress("shr_dedx_Y", &shr_dedx_Y, &b_shr_dedx_Y);
   fChain->SetBranchAddress("shr_dedx_V", &shr_dedx_V, &b_shr_dedx_V);
   fChain->SetBranchAddress("shr_dedx_U", &shr_dedx_U, &b_shr_dedx_U);
   fChain->SetBranchAddress("shr_dedx_Y_cali", &shr_dedx_Y_cali, &b_shr_dedx_Y_cali);
   fChain->SetBranchAddress("shr_dedx_V_cali", &shr_dedx_V_cali, &b_shr_dedx_V_cali);
   fChain->SetBranchAddress("shr_dedx_U_cali", &shr_dedx_U_cali, &b_shr_dedx_U_cali);
   fChain->SetBranchAddress("shr_tkfit_dedx_Y", &shr_tkfit_dedx_Y, &b_shr_tkfit_dedx_Y);
   fChain->SetBranchAddress("shr_tkfit_dedx_V", &shr_tkfit_dedx_V, &b_shr_tkfit_dedx_V);
   fChain->SetBranchAddress("shr_tkfit_dedx_U", &shr_tkfit_dedx_U, &b_shr_tkfit_dedx_U);
   fChain->SetBranchAddress("shr_tkfit_dedx_max", &shr_tkfit_dedx_max, &b_shr_tkfit_dedx_max);
   fChain->SetBranchAddress("shr_tkfit_nhits_Y", &shr_tkfit_nhits_Y, &b_shr_tkfit_nhits_Y);
   fChain->SetBranchAddress("shr_tkfit_nhits_V", &shr_tkfit_nhits_V, &b_shr_tkfit_nhits_V);
   fChain->SetBranchAddress("shr_tkfit_nhits_U", &shr_tkfit_nhits_U, &b_shr_tkfit_nhits_U);
   fChain->SetBranchAddress("shr_llrpid_dedx_Y", &shr_llrpid_dedx_Y, &b_shr_llrpid_dedx_Y);
   fChain->SetBranchAddress("shr_llrpid_dedx_V", &shr_llrpid_dedx_V, &b_shr_llrpid_dedx_V);
   fChain->SetBranchAddress("shr_llrpid_dedx_U", &shr_llrpid_dedx_U, &b_shr_llrpid_dedx_U);
   fChain->SetBranchAddress("shr_llrpid_dedx", &shr_llrpid_dedx, &b_shr_llrpid_dedx);
   fChain->SetBranchAddress("shr_tkfit_dedx_Y_alt", &shr_tkfit_dedx_Y_alt, &b_shr_tkfit_dedx_Y_alt);
   fChain->SetBranchAddress("shr_tkfit_dedx_V_alt", &shr_tkfit_dedx_V_alt, &b_shr_tkfit_dedx_V_alt);
   fChain->SetBranchAddress("shr_tkfit_dedx_U_alt", &shr_tkfit_dedx_U_alt, &b_shr_tkfit_dedx_U_alt);
   fChain->SetBranchAddress("shr_tkfit_nhits_Y_alt", &shr_tkfit_nhits_Y_alt, &b_shr_tkfit_nhits_Y_alt);
   fChain->SetBranchAddress("shr_tkfit_nhits_V_alt", &shr_tkfit_nhits_V_alt, &b_shr_tkfit_nhits_V_alt);
   fChain->SetBranchAddress("shr_tkfit_nhits_U_alt", &shr_tkfit_nhits_U_alt, &b_shr_tkfit_nhits_U_alt);
   fChain->SetBranchAddress("trkfit", &trkfit, &b__trkfit);
   fChain->SetBranchAddress("shr_tkfit_npoints", &shr_tkfit_npoints, &b_shr_tkfit_npoints);
   fChain->SetBranchAddress("shr_tkfit_npointsvalid", &shr_tkfit_npointsvalid, &b_shr_tkfit_npointsvalid);
   fChain->SetBranchAddress("shr_trkfitmedangle", &shr_trkfitmedangle, &b_shr_trkfitmedangle);
   fChain->SetBranchAddress("shrmoliereavg", &shrmoliereavg, &b_shrmoliereavg);
   fChain->SetBranchAddress("shrmoliererms", &shrmoliererms, &b_shrmoliererms);
   fChain->SetBranchAddress("shr1shr2moliereavg", &shr1shr2moliereavg, &b_shr1shr2moliereavg);
   fChain->SetBranchAddress("shr1shr2moliererms", &shr1shr2moliererms, &b_shr1shr2moliererms);
   fChain->SetBranchAddress("shr1trk1moliereavg", &shr1trk1moliereavg, &b_shr1trk1moliereavg);
   fChain->SetBranchAddress("shr1trk1moliererms", &shr1trk1moliererms, &b_shr1trk1moliererms);
   fChain->SetBranchAddress("shr1trk2moliereavg", &shr1trk2moliereavg, &b_shr1trk2moliereavg);
   fChain->SetBranchAddress("shr1trk2moliererms", &shr1trk2moliererms, &b_shr1trk2moliererms);
   fChain->SetBranchAddress("ismerged", &ismerged, &b_ismerged);
   fChain->SetBranchAddress("merge_bestdot", &merge_bestdot, &b_merge_bestdot);
   fChain->SetBranchAddress("merge_bestdist", &merge_bestdist, &b_merge_bestdist);
   fChain->SetBranchAddress("merge_vtx_x", &merge_vtx_x, &b_merge_vtx_x);
   fChain->SetBranchAddress("merge_vtx_y", &merge_vtx_y, &b_merge_vtx_y);
   fChain->SetBranchAddress("merge_vtx_z", &merge_vtx_z, &b_merge_vtx_z);
   fChain->SetBranchAddress("merge_tk_ipfp", &merge_tk_ipfp, &b_merge_tk_ipfp);
   fChain->SetBranchAddress("shr_tkfit_2cm_dedx_Y", &shr_tkfit_2cm_dedx_Y, &b_shr_tkfit_2cm_dedx_Y);
   fChain->SetBranchAddress("shr_tkfit_2cm_dedx_V", &shr_tkfit_2cm_dedx_V, &b_shr_tkfit_2cm_dedx_V);
   fChain->SetBranchAddress("shr_tkfit_2cm_dedx_U", &shr_tkfit_2cm_dedx_U, &b_shr_tkfit_2cm_dedx_U);
   fChain->SetBranchAddress("shr_tkfit_2cm_nhits_Y", &shr_tkfit_2cm_nhits_Y, &b_shr_tkfit_2cm_nhits_Y);
   fChain->SetBranchAddress("shr_tkfit_2cm_nhits_V", &shr_tkfit_2cm_nhits_V, &b_shr_tkfit_2cm_nhits_V);
   fChain->SetBranchAddress("shr_tkfit_2cm_nhits_U", &shr_tkfit_2cm_nhits_U, &b_shr_tkfit_2cm_nhits_U);
   fChain->SetBranchAddress("shr_tkfit_gap05_dedx_Y", &shr_tkfit_gap05_dedx_Y, &b_shr_tkfit_gap05_dedx_Y);
   fChain->SetBranchAddress("shr_tkfit_gap05_dedx_V", &shr_tkfit_gap05_dedx_V, &b_shr_tkfit_gap05_dedx_V);
   fChain->SetBranchAddress("shr_tkfit_gap05_dedx_U", &shr_tkfit_gap05_dedx_U, &b_shr_tkfit_gap05_dedx_U);
   fChain->SetBranchAddress("shr_tkfit_gap05_nhits_Y", &shr_tkfit_gap05_nhits_Y, &b_shr_tkfit_gap05_nhits_Y);
   fChain->SetBranchAddress("shr_tkfit_gap05_nhits_V", &shr_tkfit_gap05_nhits_V, &b_shr_tkfit_gap05_nhits_V);
   fChain->SetBranchAddress("shr_tkfit_gap05_nhits_U", &shr_tkfit_gap05_nhits_U, &b_shr_tkfit_gap05_nhits_U);
   fChain->SetBranchAddress("shr_tkfit_gap10_dedx_Y", &shr_tkfit_gap10_dedx_Y, &b_shr_tkfit_gap10_dedx_Y);
   fChain->SetBranchAddress("shr_tkfit_gap10_dedx_V", &shr_tkfit_gap10_dedx_V, &b_shr_tkfit_gap10_dedx_V);
   fChain->SetBranchAddress("shr_tkfit_gap10_dedx_U", &shr_tkfit_gap10_dedx_U, &b_shr_tkfit_gap10_dedx_U);
   fChain->SetBranchAddress("shr_tkfit_gap10_nhits_Y", &shr_tkfit_gap10_nhits_Y, &b_shr_tkfit_gap10_nhits_Y);
   fChain->SetBranchAddress("shr_tkfit_gap10_nhits_V", &shr_tkfit_gap10_nhits_V, &b_shr_tkfit_gap10_nhits_V);
   fChain->SetBranchAddress("shr_tkfit_gap10_nhits_U", &shr_tkfit_gap10_nhits_U, &b_shr_tkfit_gap10_nhits_U);
   fChain->SetBranchAddress("shr_chipr", &shr_chipr, &b_shr_chipr);
   fChain->SetBranchAddress("shr_chimu", &shr_chimu, &b_shr_chimu);
   fChain->SetBranchAddress("shr_bragg_p", &shr_bragg_p, &b_shr_bragg_p);
   fChain->SetBranchAddress("shr_bragg_mu", &shr_bragg_mu, &b_shr_bragg_mu);
   fChain->SetBranchAddress("shr_bragg_mip", &shr_bragg_mip, &b_shr_bragg_mip);
   fChain->SetBranchAddress("shr_bragg_kaon", &shr_bragg_kaon, &b_shr_bragg_kaon);
   fChain->SetBranchAddress("shr_bragg_pion", &shr_bragg_pion, &b_shr_bragg_pion);
   fChain->SetBranchAddress("tksh_distance", &tksh_distance, &b_tksh_distance);
   fChain->SetBranchAddress("tksh_angle", &tksh_angle, &b_tksh_angle);
   fChain->SetBranchAddress("shr_distance", &shr_distance, &b_shr_distance);
   fChain->SetBranchAddress("shr_score", &shr_score, &b_shr_score);
   fChain->SetBranchAddress("shr_bkt_pdg", &shr_bkt_pdg, &b_shr_bkt_pdg);
   fChain->SetBranchAddress("shr_bkt_purity", &shr_bkt_purity, &b_shr_bkt_purity);
   fChain->SetBranchAddress("shr_bkt_completeness", &shr_bkt_completeness, &b_shr_bkt_completeness);
   fChain->SetBranchAddress("shr_bkt_E", &shr_bkt_E, &b_shr_bkt_E);
   fChain->SetBranchAddress("trk_len", &trk_len, &b_trk_len);
   fChain->SetBranchAddress("trk_theta", &trk_theta, &b_trk_theta);
   fChain->SetBranchAddress("trk_phi", &trk_phi, &b_trk_phi);
   fChain->SetBranchAddress("trk_energy", &trk_energy, &b_trk_energy);
   fChain->SetBranchAddress("trk_energy_muon", &trk_energy_muon, &b_trk_energy_muon);
   fChain->SetBranchAddress("trk_energy_muon_mcs", &trk_energy_muon_mcs, &b_trk_energy_muon_mcs);
   fChain->SetBranchAddress("trk_energy_tot", &trk_energy_tot, &b_trk_energy_tot);
   fChain->SetBranchAddress("trk_energy_muon_tot", &trk_energy_muon_tot, &b_trk_energy_muon_tot);
   fChain->SetBranchAddress("trk_distance", &trk_distance, &b_trk_distance);
   fChain->SetBranchAddress("trk_score", &trk_score, &b_trk_score);
   fChain->SetBranchAddress("trk_bkt_pdg", &trk_bkt_pdg, &b_trk_bkt_pdg);
   fChain->SetBranchAddress("trk_bkt_purity", &trk_bkt_purity, &b_trk_bkt_purity);
   fChain->SetBranchAddress("trk_bkt_completeness", &trk_bkt_completeness, &b_trk_bkt_completeness);
   fChain->SetBranchAddress("trk_bkt_E", &trk_bkt_E, &b_trk_bkt_E);
   fChain->SetBranchAddress("trk_chipr_best", &trk_chipr_best, &b_trk_chipr_best);
   fChain->SetBranchAddress("trk_chipr_worst", &trk_chipr_worst, &b_trk_chipr_worst);
   fChain->SetBranchAddress("trk_chimu_best", &trk_chimu_best, &b_trk_chimu_best);
   fChain->SetBranchAddress("trk_chimu_worst", &trk_chimu_worst, &b_trk_chimu_worst);
   fChain->SetBranchAddress("trk_chipr", &trk_chipr, &b_trk_chipr);
   fChain->SetBranchAddress("trk_chimu", &trk_chimu, &b_trk_chimu);
   fChain->SetBranchAddress("trk_pida", &trk_pida, &b_trk_pida);
   fChain->SetBranchAddress("trk_bragg_p", &trk_bragg_p, &b_trk_bragg_p);
   fChain->SetBranchAddress("trk_bragg_mu", &trk_bragg_mu, &b_trk_bragg_mu);
   fChain->SetBranchAddress("trk_bragg_mip", &trk_bragg_mip, &b_trk_bragg_mip);
   fChain->SetBranchAddress("trk_bragg_kaon", &trk_bragg_kaon, &b_trk_bragg_kaon);
   fChain->SetBranchAddress("trk_bragg_pion", &trk_bragg_pion, &b_trk_bragg_pion);
   fChain->SetBranchAddress("trk_hits_max", &trk_hits_max, &b_trk_hits_max);
   fChain->SetBranchAddress("shr_hits_max", &shr_hits_max, &b_shr_hits_max);
   fChain->SetBranchAddress("all_shr_hits", &all_shr_hits, &b_all_shr_hits);
   fChain->SetBranchAddress("all_trk_hits", &all_trk_hits, &b_all_trk_hits);
   fChain->SetBranchAddress("all_shr_energies", &all_shr_energies, &b_all_shr_energies);
   fChain->SetBranchAddress("all_trk_energies", &all_trk_energies, &b_all_trk_energies);
//    fChain->SetBranchAddress("shr_hits_max", &shr_hits_max, &b_shr_hits_max);
   fChain->SetBranchAddress("trk_hits_2nd", &trk_hits_2nd, &b_trk_hits_2nd);
   fChain->SetBranchAddress("shr_hits_2nd", &shr_hits_2nd, &b_shr_hits_2nd);
   fChain->SetBranchAddress("trk_hits_3rd", &trk_hits_3rd, &b_trk_hits_3rd);
   fChain->SetBranchAddress("shr_hits_3rd", &shr_hits_3rd, &b_shr_hits_3rd);
   fChain->SetBranchAddress("trkshrhitdist0", &trkshrhitdist0, &b_trkshrhitdist0);
   fChain->SetBranchAddress("trkshrhitdist1", &trkshrhitdist1, &b_trkshrhitdist1);
   fChain->SetBranchAddress("trkshrhitdist2", &trkshrhitdist2, &b_trkshrhitdist2);
   fChain->SetBranchAddress("trk2shrhitdist0", &trk2shrhitdist0, &b_trk2shrhitdist0);
   fChain->SetBranchAddress("trk2shrhitdist1", &trk2shrhitdist1, &b_trk2shrhitdist1);
   fChain->SetBranchAddress("trk2shrhitdist2", &trk2shrhitdist2, &b_trk2shrhitdist2);
   fChain->SetBranchAddress("trk1trk2hitdist0", &trk1trk2hitdist0, &b_trk1trk2hitdist0);
   fChain->SetBranchAddress("trk1trk2hitdist1", &trk1trk2hitdist1, &b_trk1trk2hitdist1);
   fChain->SetBranchAddress("trk1trk2hitdist2", &trk1trk2hitdist2, &b_trk1trk2hitdist2);
   fChain->SetBranchAddress("total_hits_y", &total_hits_y, &b_total_hits_y);
   fChain->SetBranchAddress("extra_energy_y", &extra_energy_y, &b_extra_energy_y);
   fChain->SetBranchAddress("trk_energy_hits_tot", &trk_energy_hits_tot, &b_trk_energy_hits_tot);
   fChain->SetBranchAddress("subcluster", &subcluster, &b_subcluster);
   fChain->SetBranchAddress("shrsubclusters0", &shrsubclusters0, &b_shrsubclusters0);
   fChain->SetBranchAddress("shrsubclusters1", &shrsubclusters1, &b_shrsubclusters1);
   fChain->SetBranchAddress("shrsubclusters2", &shrsubclusters2, &b_shrsubclusters2);
   fChain->SetBranchAddress("shrclusfrac0", &shrclusfrac0, &b_shrclusfrac0);
   fChain->SetBranchAddress("shrclusfrac1", &shrclusfrac1, &b_shrclusfrac1);
   fChain->SetBranchAddress("shrclusfrac2", &shrclusfrac2, &b_shrclusfrac2);
   fChain->SetBranchAddress("shrclusdir0", &shrclusdir0, &b_shrclusdir0);
   fChain->SetBranchAddress("shrclusdir1", &shrclusdir1, &b_shrclusdir1);
   fChain->SetBranchAddress("shrclusdir2", &shrclusdir2, &b_shrclusdir2);
   fChain->SetBranchAddress("shr_hits_tot", &shr_hits_tot, &b_shr_hits_tot);
   fChain->SetBranchAddress("shr_hits_y_tot", &shr_hits_y_tot, &b_shr_hits_y_tot);
   fChain->SetBranchAddress("shr_hits_u_tot", &shr_hits_u_tot, &b_shr_hits_u_tot);
   fChain->SetBranchAddress("shr_hits_v_tot", &shr_hits_v_tot, &b_shr_hits_v_tot);
   fChain->SetBranchAddress("trk_hits_tot", &trk_hits_tot, &b_trk_hits_tot);
   fChain->SetBranchAddress("trk_hits_y_tot", &trk_hits_y_tot, &b_trk_hits_y_tot);
   fChain->SetBranchAddress("trk_hits_u_tot", &trk_hits_u_tot, &b_trk_hits_u_tot);
   fChain->SetBranchAddress("trk_hits_v_tot", &trk_hits_v_tot, &b_trk_hits_v_tot);
   fChain->SetBranchAddress("_elecclusters_U_charge", &_elecclusters_U_charge, &b_elecclusters_U_charge);
   fChain->SetBranchAddress("_elecclusters_V_charge", &_elecclusters_V_charge, &b_elecclusters_V_charge);
   fChain->SetBranchAddress("_elecclusters_Y_charge", &_elecclusters_Y_charge, &b_elecclusters_Y_charge);
   fChain->SetBranchAddress("_elecclusters_U_N", &_elecclusters_U_N, &b_elecclusters_U_N);
   fChain->SetBranchAddress("_elecclusters_V_N", &_elecclusters_V_N, &b_elecclusters_V_N);
   fChain->SetBranchAddress("_elecclusters_Y_N", &_elecclusters_Y_N, &b_elecclusters_Y_N);
   fChain->SetBranchAddress("n_tracks_contained", &n_tracks_contained, &b_n_tracks_contained);
   fChain->SetBranchAddress("n_showers_contained", &n_showers_contained, &b_n_showers_contained);
   fChain->SetBranchAddress("matched_E", &matched_E, &b_matched_E);
   fChain->SetBranchAddress("hits_ratio", &hits_ratio, &b_hits_ratio);
   fChain->SetBranchAddress("contained_fraction", &contained_fraction, &b_contained_fraction);
   fChain->SetBranchAddress("sps_contained_fraction", &sps_contained_fraction, &b_sps_contained_fraction);
   fChain->SetBranchAddress("pt", &pt, &b_pt);
   fChain->SetBranchAddress("p", &p, &b_p);
   fChain->SetBranchAddress("pt_assume_muon", &pt_assume_muon, &b_pt_assume_muon);
   fChain->SetBranchAddress("p_assume_muon", &p_assume_muon, &b_p_assume_muon);
   fChain->SetBranchAddress("reco_e", &reco_e, &b_reco_e);
   fChain->SetBranchAddress("dvtx", &dvtx, &b_dvtx);
   fChain->SetBranchAddress("dtrk", &dtrk, &b_dtrk);
   fChain->SetBranchAddress("contained_sps_ratio", &contained_sps_ratio, &b_contained_sps_ratio);
   fChain->SetBranchAddress("dtrk_x_boundary", &dtrk_x_boundary, &b_dtrk_x_boundary);
   fChain->SetBranchAddress("dtrk_y_boundary", &dtrk_y_boundary, &b_dtrk_y_boundary);
   fChain->SetBranchAddress("dtrk_z_boundary", &dtrk_z_boundary, &b_dtrk_z_boundary);
   fChain->SetBranchAddress("dshr_x_boundary", &dshr_x_boundary, &b_dshr_x_boundary);
   fChain->SetBranchAddress("dshr_y_boundary", &dshr_y_boundary, &b_dshr_y_boundary);
   fChain->SetBranchAddress("dshr_z_boundary", &dshr_z_boundary, &b_dshr_z_boundary);
   fChain->SetBranchAddress("dvtx_x_boundary", &dvtx_x_boundary, &b_dvtx_x_boundary);
   fChain->SetBranchAddress("dvtx_y_boundary", &dvtx_y_boundary, &b_dvtx_y_boundary);
   fChain->SetBranchAddress("dvtx_z_boundary", &dvtx_z_boundary, &b_dvtx_z_boundary);
   fChain->SetBranchAddress("dtrk_boundary", &dtrk_boundary, &b_dtrk_boundary);
   fChain->SetBranchAddress("dvtx_boundary", &dvtx_boundary, &b_dvtx_boundary);
   fChain->SetBranchAddress("dshr_boundary", &dshr_boundary, &b_dshr_boundary);
   fChain->SetBranchAddress("dmc_boundary", &dmc_boundary, &b_dmc_boundary);
   fChain->SetBranchAddress("CosmicIP", &CosmicIP, &b_CosmicIP);
   fChain->SetBranchAddress("CosmicIPAll3D", &CosmicIPAll3D, &b_CosmicIPAll3D);
   fChain->SetBranchAddress("CosmicDirAll3D", &CosmicDirAll3D, &b_CosmicDirAll3D);
   fChain->SetBranchAddress("CosmicIPAll2DEnds", &CosmicIPAll2DEnds, &b_CosmicIPAll2DEnds);
   fChain->SetBranchAddress("CosmicDirAll2DEnds", &CosmicDirAll2DEnds, &b_CosmicDirAll2DEnds);
   fChain->SetBranchAddress("CosmicIPAll2DOvlp", &CosmicIPAll2DOvlp, &b_CosmicIPAll2DOvlp);
   fChain->SetBranchAddress("CosmicDirAll2DOvlp", &CosmicDirAll2DOvlp, &b_CosmicDirAll2DOvlp);
   fChain->SetBranchAddress("leeweight", &leeweight, &b_leeweight);
   fChain->SetBranchAddress("true_pt", &true_pt, &b_true_pt);
   fChain->SetBranchAddress("true_pt_visible", &true_pt_visible, &b_true_pt_visible);
   fChain->SetBranchAddress("true_p", &true_p, &b_true_p);
   fChain->SetBranchAddress("true_p_visible", &true_p_visible, &b_true_p_visible);
   fChain->SetBranchAddress("true_e_visible", &true_e_visible, &b_true_e_visible);
   fChain->SetBranchAddress("_opfilter_pe_beam", &_opfilter_pe_beam, &b_opfilter_pe_beam);
   fChain->SetBranchAddress("_opfilter_pe_veto", &_opfilter_pe_veto, &b_opfilter_pe_veto);
   fChain->SetBranchAddress("nu_pdg", &nu_pdg, &b_nu_pdg);
   fChain->SetBranchAddress("ccnc", &ccnc, &b_ccnc);
   fChain->SetBranchAddress("nu_parent_pdg", &nu_parent_pdg, &b_nu_parent_pdg);
   fChain->SetBranchAddress("nu_hadron_pdg", &nu_hadron_pdg, &b_nu_hadron_pdg);
   fChain->SetBranchAddress("nu_decay_mode", &nu_decay_mode, &b_nu_decay_mode);
   fChain->SetBranchAddress("par_decay_vx", &par_decay_vx, &b_par_decay_vx);
   fChain->SetBranchAddress("par_decay_vy", &par_decay_vy, &b_par_decay_vy);
   fChain->SetBranchAddress("par_decay_vz", &par_decay_vz, &b_par_decay_vz);
   fChain->SetBranchAddress("par_decay_px", &par_decay_px, &b_par_decay_px);
   fChain->SetBranchAddress("par_decay_py", &par_decay_py, &b_par_decay_py);
   fChain->SetBranchAddress("par_decay_pz", &par_decay_pz, &b_par_decay_pz);
   fChain->SetBranchAddress("baseline", &baseline, &b_baseline);
   fChain->SetBranchAddress("true_nu_px", &true_nu_px, &b_true_nu_px);
   fChain->SetBranchAddress("true_nu_py", &true_nu_py, &b_true_nu_py);
   fChain->SetBranchAddress("true_nu_pz", &true_nu_pz, &b_true_nu_pz);
   fChain->SetBranchAddress("swtrig_pre", &swtrig_pre, &b_swtrig_pre);
   fChain->SetBranchAddress("swtrig_post", &swtrig_post, &b_swtrig_post);
   fChain->SetBranchAddress("swtrig_pre_ext", &swtrig_pre_ext, &b_swtrig_pre_ext);
   fChain->SetBranchAddress("swtrig_post_ext", &swtrig_post_ext, &b_swtrig_post_ext);
   fChain->SetBranchAddress("interaction", &interaction, &b_interaction);
   fChain->SetBranchAddress("nu_e", &nu_e, &b_nu_e);
   fChain->SetBranchAddress("nu_l", &nu_l, &b_nu_l);
   fChain->SetBranchAddress("nu_pt", &nu_pt, &b_nu_pt);
   fChain->SetBranchAddress("theta", &theta, &b_theta);
   fChain->SetBranchAddress("isVtxInFiducial", &isVtxInFiducial, &b_isVtxInFiducial);
   fChain->SetBranchAddress("truthFiducial", &truthFiducial, &b_truthFiducial);
   fChain->SetBranchAddress("true_nu_vtx_t", &true_nu_vtx_t, &b_true_nu_vtx_t);
   fChain->SetBranchAddress("true_nu_vtx_x", &true_nu_vtx_x, &b_true_nu_vtx_x);
   fChain->SetBranchAddress("true_nu_vtx_y", &true_nu_vtx_y, &b_true_nu_vtx_y);
   fChain->SetBranchAddress("true_nu_vtx_z", &true_nu_vtx_z, &b_true_nu_vtx_z);
   fChain->SetBranchAddress("true_nu_vtx_sce_x", &true_nu_vtx_sce_x, &b_true_nu_vtx_sce_x);
   fChain->SetBranchAddress("true_nu_vtx_sce_y", &true_nu_vtx_sce_y, &b_true_nu_vtx_sce_y);
   fChain->SetBranchAddress("true_nu_vtx_sce_z", &true_nu_vtx_sce_z, &b_true_nu_vtx_sce_z);
   fChain->SetBranchAddress("reco_nu_vtx_x", &reco_nu_vtx_x, &b_reco_nu_vtx_x);
   fChain->SetBranchAddress("reco_nu_vtx_y", &reco_nu_vtx_y, &b_reco_nu_vtx_y);
   fChain->SetBranchAddress("reco_nu_vtx_z", &reco_nu_vtx_z, &b_reco_nu_vtx_z);
   fChain->SetBranchAddress("reco_nu_vtx_sce_x", &reco_nu_vtx_sce_x, &b_reco_nu_vtx_sce_x);
   fChain->SetBranchAddress("reco_nu_vtx_sce_y", &reco_nu_vtx_sce_y, &b_reco_nu_vtx_sce_y);
   fChain->SetBranchAddress("reco_nu_vtx_sce_z", &reco_nu_vtx_sce_z, &b_reco_nu_vtx_sce_z);
   fChain->SetBranchAddress("nmuon", &nmuon, &b_nmuon);
   fChain->SetBranchAddress("muon_e", &muon_e, &b_muon_e);
   fChain->SetBranchAddress("muon_c", &muon_c, &b_muon_c);
   fChain->SetBranchAddress("muon_p", &muon_p, &b_muon_p);
   fChain->SetBranchAddress("nelec", &nelec, &b_nelec);
   fChain->SetBranchAddress("elec_e", &elec_e, &b_elec_e);
   fChain->SetBranchAddress("elec_c", &elec_c, &b_elec_c);
   fChain->SetBranchAddress("elec_p", &elec_p, &b_elec_p);
   fChain->SetBranchAddress("elec_vx", &elec_vx, &b_elec_vx);
   fChain->SetBranchAddress("elec_vy", &elec_vy, &b_elec_vy);
   fChain->SetBranchAddress("elec_vz", &elec_vz, &b_elec_vz);
   fChain->SetBranchAddress("elec_px", &elec_px, &b_elec_px);
   fChain->SetBranchAddress("elec_py", &elec_py, &b_elec_py);
   fChain->SetBranchAddress("elec_pz", &elec_pz, &b_elec_pz);
   fChain->SetBranchAddress("npi0", &npi0, &b_npi0);
   fChain->SetBranchAddress("pi0_e", &pi0_e, &b_pi0_e);
   fChain->SetBranchAddress("pi0_c", &pi0_c, &b_pi0_c);
   fChain->SetBranchAddress("pi0_p", &pi0_p, &b_pi0_p);
   fChain->SetBranchAddress("nneutron", &nneutron, &b_nneutron);
   fChain->SetBranchAddress("nproton", &nproton, &b_nproton);
   fChain->SetBranchAddress("proton_e", &proton_e, &b_proton_e);
   fChain->SetBranchAddress("proton_c", &proton_c, &b_proton_c);
   fChain->SetBranchAddress("proton_p", &proton_p, &b_proton_p);
   fChain->SetBranchAddress("npion", &npion, &b_npion);
   fChain->SetBranchAddress("pion_e", &pion_e, &b_pion_e);
   fChain->SetBranchAddress("pion_c", &pion_c, &b_pion_c);
   fChain->SetBranchAddress("pion_p", &pion_p, &b_pion_p);
   fChain->SetBranchAddress("neta", &neta, &b_neta);
   fChain->SetBranchAddress("eta_e", &eta_e, &b_eta_e);
   fChain->SetBranchAddress("nslice", &nslice, &b_nslice);
   fChain->SetBranchAddress("crtveto", &crtveto, &b_crtveto);
   fChain->SetBranchAddress("crthitpe", &crthitpe, &b_crthitpe);
   fChain->SetBranchAddress("pfp_slice_idx", &pfp_slice_idx, &b_pfp_slice_idx);
   fChain->SetBranchAddress("category", &category, &b_category);
   fChain->SetBranchAddress("backtracked_pdg", &backtracked_pdg, &b_backtracked_pdg);
   fChain->SetBranchAddress("backtracked_e", &backtracked_e, &b_backtracked_e);
   fChain->SetBranchAddress("backtracked_tid", &backtracked_tid, &b_backtracked_tid);
   fChain->SetBranchAddress("backtracked_purity", &backtracked_purity, &b_backtracked_purity);
   fChain->SetBranchAddress("backtracked_completeness", &backtracked_completeness, &b_backtracked_completeness);
   fChain->SetBranchAddress("backtracked_overlay_purity", &backtracked_overlay_purity, &b_backtracked_overlay_purity);
   fChain->SetBranchAddress("backtracked_px", &backtracked_px, &b_backtracked_px);
   fChain->SetBranchAddress("backtracked_py", &backtracked_py, &b_backtracked_py);
   fChain->SetBranchAddress("backtracked_pz", &backtracked_pz, &b_backtracked_pz);
   fChain->SetBranchAddress("backtracked_start_x", &backtracked_start_x, &b_backtracked_start_x);
   fChain->SetBranchAddress("backtracked_start_y", &backtracked_start_y, &b_backtracked_start_y);
   fChain->SetBranchAddress("backtracked_start_z", &backtracked_start_z, &b_backtracked_start_z);
   fChain->SetBranchAddress("backtracked_start_t", &backtracked_start_t, &b_backtracked_start_t);
   fChain->SetBranchAddress("backtracked_start_U", &backtracked_start_U, &b_backtracked_start_U);
   fChain->SetBranchAddress("backtracked_start_V", &backtracked_start_V, &b_backtracked_start_V);
   fChain->SetBranchAddress("backtracked_start_Y", &backtracked_start_Y, &b_backtracked_start_Y);
   fChain->SetBranchAddress("backtracked_sce_start_x", &backtracked_sce_start_x, &b_backtracked_sce_start_x);
   fChain->SetBranchAddress("backtracked_sce_start_y", &backtracked_sce_start_y, &b_backtracked_sce_start_y);
   fChain->SetBranchAddress("backtracked_sce_start_z", &backtracked_sce_start_z, &b_backtracked_sce_start_z);
   fChain->SetBranchAddress("backtracked_sce_start_U", &backtracked_sce_start_U, &b_backtracked_sce_start_U);
   fChain->SetBranchAddress("backtracked_sce_start_V", &backtracked_sce_start_V, &b_backtracked_sce_start_V);
   fChain->SetBranchAddress("backtracked_sce_start_Y", &backtracked_sce_start_Y, &b_backtracked_sce_start_Y);
   fChain->SetBranchAddress("lep_e", &lep_e, &b_lep_e);
   fChain->SetBranchAddress("pass", &pass, &b_pass);
   fChain->SetBranchAddress("swtrig", &swtrig, &b_swtrig);
   fChain->SetBranchAddress("evnhits", &evnhits, &b_evnhits);
   fChain->SetBranchAddress("slpdg", &slpdg, &b_slpdg);
   fChain->SetBranchAddress("slnhits", &slnhits, &b_slnhits);
   fChain->SetBranchAddress("n_pfps", &n_pfps, &b_n_pfps);
   fChain->SetBranchAddress("n_tracks", &n_tracks, &b_n_tracks);
   fChain->SetBranchAddress("n_showers", &n_showers, &b_n_showers);
   fChain->SetBranchAddress("pfp_generation_v", &pfp_generation_v, &b_pfp_generation_v);
   fChain->SetBranchAddress("pfp_trk_daughters_v", &pfp_trk_daughters_v, &b_pfp_trk_daughters_v);
   fChain->SetBranchAddress("pfp_shr_daughters_v", &pfp_shr_daughters_v, &b_pfp_shr_daughters_v);
   fChain->SetBranchAddress("pfp_n_descendents_v", &pfp_n_descendents_v, &b_pfp_n_descendents_v);
   fChain->SetBranchAddress("pfp_vtx_x_v", &pfp_vtx_x_v, &b_pfp_vtx_x_v);
   fChain->SetBranchAddress("pfp_vtx_y_v", &pfp_vtx_y_v, &b_pfp_vtx_y_v);
   fChain->SetBranchAddress("pfp_vtx_z_v", &pfp_vtx_z_v, &b_pfp_vtx_z_v);
   fChain->SetBranchAddress("trk_score_v", &trk_score_v, &b_trk_score_v);
   fChain->SetBranchAddress("pfpdg", &pfpdg, &b_pfpdg);
   fChain->SetBranchAddress("pfnhits", &pfnhits, &b_pfnhits);
   fChain->SetBranchAddress("pfnplanehits_U", &pfnplanehits_U, &b_pfnplanehits_U);
   fChain->SetBranchAddress("pfnplanehits_V", &pfnplanehits_V, &b_pfnplanehits_V);
   fChain->SetBranchAddress("pfnplanehits_Y", &pfnplanehits_Y, &b_pfnplanehits_Y);
   fChain->SetBranchAddress("pfpplanesubclusters_U", &pfpplanesubclusters_U, &b_pfpplanesubclusters_U);
   fChain->SetBranchAddress("pfpplanesubclusters_V", &pfpplanesubclusters_V, &b_pfpplanesubclusters_V);
   fChain->SetBranchAddress("pfpplanesubclusters_Y", &pfpplanesubclusters_Y, &b_pfpplanesubclusters_Y);
   fChain->SetBranchAddress("pfpplanesubhitfracmax_U", &pfpplanesubhitfracmax_U, &b_pfpplanesubhitfracmax_U);
   fChain->SetBranchAddress("pfpplanesubhitfracmax_V", &pfpplanesubhitfracmax_V, &b_pfpplanesubhitfracmax_V);
   fChain->SetBranchAddress("pfpplanesubhitfracmax_Y", &pfpplanesubhitfracmax_Y, &b_pfpplanesubhitfracmax_Y);
   fChain->SetBranchAddress("hits_u", &hits_u, &b_hits_u);
   fChain->SetBranchAddress("hits_v", &hits_v, &b_hits_v);
   fChain->SetBranchAddress("hits_y", &hits_y, &b_hits_y);
   fChain->SetBranchAddress("slice_id", &slice_id, &b_slice_id);
   fChain->SetBranchAddress("slice_topo_score_v", &slice_topo_score_v, &b_slice_topo_score_v);
   fChain->SetBranchAddress("topological_score", &topological_score, &b_topological_score);
   fChain->SetBranchAddress("slclustfrac", &slclustfrac, &b_slclustfrac);
   fChain->SetBranchAddress("mc_pdg", &mc_pdg, &b_mc_pdg);
   fChain->SetBranchAddress("mc_E", &mc_E, &b_mc_E);
   fChain->SetBranchAddress("mc_n_elastic", &mc_n_elastic, &b_mc_n_elastic);
   fChain->SetBranchAddress("mc_n_inelastic", &mc_n_inelastic, &b_mc_n_inelastic);
   fChain->SetBranchAddress("mc_vx", &mc_vx, &b_mc_vx);
   fChain->SetBranchAddress("mc_vy", &mc_vy, &b_mc_vy);
   fChain->SetBranchAddress("mc_vz", &mc_vz, &b_mc_vz);
   fChain->SetBranchAddress("mc_endx", &mc_endx, &b_mc_endx);
   fChain->SetBranchAddress("mc_endy", &mc_endy, &b_mc_endy);
   fChain->SetBranchAddress("mc_endz", &mc_endz, &b_mc_endz);
   fChain->SetBranchAddress("mc_px", &mc_px, &b_mc_px);
   fChain->SetBranchAddress("mc_py", &mc_py, &b_mc_py);
   fChain->SetBranchAddress("mc_pz", &mc_pz, &b_mc_pz);
   fChain->SetBranchAddress("mc_end_p", &mc_end_p, &b_mc_end_p);
   fChain->SetBranchAddress("mc_completeness", &mc_completeness, &b_mc_completeness);
   fChain->SetBranchAddress("mc_purity", &mc_purity, &b_mc_purity);
   fChain->SetBranchAddress("endmuonprocess", &endmuonprocess, &b_endmuonprocess);
   fChain->SetBranchAddress("endmuonmichel", &endmuonmichel, &b_endmuonmichel);
   fChain->SetBranchAddress("filter_antibdt", &filter_antibdt, &b_filter_antibdt);
   fChain->SetBranchAddress("filter_ncpi0", &filter_ncpi0, &b_filter_ncpi0);
   fChain->SetBranchAddress("filter_pi0", &filter_pi0, &b_filter_pi0);
   fChain->SetBranchAddress("filter_ccinclusive", &filter_ccinclusive, &b_filter_ccinclusive);
   fChain->SetBranchAddress("weights", &weights, &b_weights);
   fChain->SetBranchAddress("weightsFlux", &weightsFlux, &b_weightsFlux);
   fChain->SetBranchAddress("weightsGenie", &weightsGenie, &b_weightsGenie);
   fChain->SetBranchAddress("weightsReint", &weightsReint, &b_weightsReint);
   fChain->SetBranchAddress("weightSpline", &weightSpline, &b_weightSpline);
   fChain->SetBranchAddress("weightTune", &weightTune, &b_weightTune);
   fChain->SetBranchAddress("weightSplineTimesTune", &weightSplineTimesTune, &b_weightSplineTimesTune);
   fChain->SetBranchAddress("ppfx_cv", &ppfx_cv, &b_ppfx_cv);
   fChain->SetBranchAddress("weightsPPFX", &weightsPPFX, &b_weightsPPFX);
   fChain->SetBranchAddress("knobRPAup", &knobRPAup, &b_knobRPAup);
   fChain->SetBranchAddress("knobRPAdn", &knobRPAdn, &b_knobRPAdn);
   fChain->SetBranchAddress("knobCCMECup", &knobCCMECup, &b_knobCCMECup);
   fChain->SetBranchAddress("knobCCMECdn", &knobCCMECdn, &b_knobCCMECdn);
   fChain->SetBranchAddress("knobAxFFCCQEup", &knobAxFFCCQEup, &b_knobAxFFCCQEup);
   fChain->SetBranchAddress("knobAxFFCCQEdn", &knobAxFFCCQEdn, &b_knobAxFFCCQEdn);
   fChain->SetBranchAddress("knobVecFFCCQEup", &knobVecFFCCQEup, &b_knobVecFFCCQEup);
   fChain->SetBranchAddress("knobVecFFCCQEdn", &knobVecFFCCQEdn, &b_knobVecFFCCQEdn);
   fChain->SetBranchAddress("knobDecayAngMECup", &knobDecayAngMECup, &b_knobDecayAngMECup);
   fChain->SetBranchAddress("knobDecayAngMECdn", &knobDecayAngMECdn, &b_knobDecayAngMECdn);
   fChain->SetBranchAddress("knobThetaDelta2Npiup", &knobThetaDelta2Npiup, &b_knobThetaDelta2Npiup);
   fChain->SetBranchAddress("knobThetaDelta2Npidn", &knobThetaDelta2Npidn, &b_knobThetaDelta2Npidn);
   fChain->SetBranchAddress("knobThetaDelta2NRadup", &knobThetaDelta2NRadup, &b_knobThetaDelta2NRadup);
   fChain->SetBranchAddress("knobThetaDelta2NRaddn", &knobThetaDelta2NRaddn, &b_knobThetaDelta2NRaddn);
   fChain->SetBranchAddress("knobNormCCCOHup", &knobNormCCCOHup, &b_knobNormCCCOHup);
   fChain->SetBranchAddress("knobNormCCCOHdn", &knobNormCCCOHdn, &b_knobNormCCCOHdn);
   fChain->SetBranchAddress("knobNormNCCOHup", &knobNormNCCOHup, &b_knobNormNCCOHup);
   fChain->SetBranchAddress("knobNormNCCOHdn", &knobNormNCCOHdn, &b_knobNormNCCOHdn);
   fChain->SetBranchAddress("knobxsr_scc_Fv3up", &knobxsr_scc_Fv3up, &b_knobxsr_scc_Fv3up);
   fChain->SetBranchAddress("knobxsr_scc_Fv3dn", &knobxsr_scc_Fv3dn, &b_knobxsr_scc_Fv3dn);
   fChain->SetBranchAddress("knobxsr_scc_Fa3up", &knobxsr_scc_Fa3up, &b_knobxsr_scc_Fa3up);
   fChain->SetBranchAddress("knobxsr_scc_Fa3dn", &knobxsr_scc_Fa3dn, &b_knobxsr_scc_Fa3dn);
   fChain->SetBranchAddress("RootinoFix", &RootinoFix, &b_RootinoFix);
   fChain->SetBranchAddress("flash_pe_flash_matching", &flash_pe_flash_matching, &b_flash_pe_flash_matching);
   fChain->SetBranchAddress("flash_pe_flash_matching_v", &flash_pe_flash_matching_v, &b_flash_pe_flash_matching_v);
   fChain->SetBranchAddress("slice_pe_flash_matching_v", &slice_pe_flash_matching_v, &b_slice_pe_flash_matching_v);
   fChain->SetBranchAddress("flash_time_flash_matching", &flash_time_flash_matching, &b_flash_time_flash_matching);
   fChain->SetBranchAddress("flash_y_flash_matching", &flash_y_flash_matching, &b_flash_y_flash_matching);
   fChain->SetBranchAddress("flash_z_flash_matching", &flash_z_flash_matching, &b_flash_z_flash_matching);
   fChain->SetBranchAddress("flash_timewidth_flash_matching", &flash_timewidth_flash_matching, &b_flash_timewidth_flash_matching);
   fChain->SetBranchAddress("flash_ywidth_flash_matching", &flash_ywidth_flash_matching, &b_flash_ywidth_flash_matching);
   fChain->SetBranchAddress("flash_zwidth_flash_matching", &flash_zwidth_flash_matching, &b_flash_zwidth_flash_matching);
   fChain->SetBranchAddress("nu_flashmatch_score", &nu_flashmatch_score, &b_nu_flashmatch_score);
   fChain->SetBranchAddress("nu_centerX", &nu_centerX, &b_nu_centerX);
   fChain->SetBranchAddress("nu_centerY", &nu_centerY, &b_nu_centerY);
   fChain->SetBranchAddress("nu_centerZ", &nu_centerZ, &b_nu_centerZ);
   fChain->SetBranchAddress("nu_totalCharge", &nu_totalCharge, &b_nu_totalCharge);
   fChain->SetBranchAddress("best_cosmic_flashmatch_score", &best_cosmic_flashmatch_score, &b_best_cosmic_flashmatch_score);
   fChain->SetBranchAddress("best_obviouscosmic_flashmatch_score", &best_obviouscosmic_flashmatch_score, &b_best_obviouscosmic_flashmatch_score);
   fChain->SetBranchAddress("cosmic_flashmatch_score_v", &cosmic_flashmatch_score_v, &b_cosmic_flashmatch_score_v);
   fChain->SetBranchAddress("cosmic_topological_score_v", &cosmic_topological_score_v, &b_cosmic_topological_score_v);
   fChain->SetBranchAddress("cosmic_centerX_v", &cosmic_centerX_v, &b_cosmic_centerX_v);
   fChain->SetBranchAddress("cosmic_centerY_v", &cosmic_centerY_v, &b_cosmic_centerY_v);
   fChain->SetBranchAddress("cosmic_centerZ_v", &cosmic_centerZ_v, &b_cosmic_centerZ_v);
   fChain->SetBranchAddress("cosmic_totalCharge_v", &cosmic_totalCharge_v, &b_cosmic_totalCharge_v);
   fChain->SetBranchAddress("cosmic_nhits_v", &cosmic_nhits_v, &b_cosmic_nhits_v);
   fChain->SetBranchAddress("cosmic_nunhits_v", &cosmic_nunhits_v, &b_cosmic_nunhits_v);
   fChain->SetBranchAddress("cosmic_isclear_v", &cosmic_isclear_v, &b_cosmic_isclear_v);
   fChain->SetBranchAddress("mcf_nu_e", &mcf_nu_e, &b_mcf_nu_e);
   fChain->SetBranchAddress("mcf_lep_e", &mcf_lep_e, &b_mcf_lep_e);
   fChain->SetBranchAddress("mcf_actvol", &mcf_actvol, &b_mcf_actvol);
   fChain->SetBranchAddress("mcf_nmm", &mcf_nmm, &b_mcf_nmm);
   fChain->SetBranchAddress("mcf_nmp", &mcf_nmp, &b_mcf_nmp);
   fChain->SetBranchAddress("mcf_nem", &mcf_nem, &b_mcf_nem);
   fChain->SetBranchAddress("mcf_nep", &mcf_nep, &b_mcf_nep);
   fChain->SetBranchAddress("mcf_np0", &mcf_np0, &b_mcf_np0);
   fChain->SetBranchAddress("mcf_npp", &mcf_npp, &b_mcf_npp);
   fChain->SetBranchAddress("mcf_npm", &mcf_npm, &b_mcf_npm);
   fChain->SetBranchAddress("mcf_mcshr_elec_etot", &mcf_mcshr_elec_etot, &b_mcf_mcshr_elec_etot);
   fChain->SetBranchAddress("mcf_pass_ccpi0", &mcf_pass_ccpi0, &b_mcf_pass_ccpi0);
   fChain->SetBranchAddress("mcf_pass_ncpi0", &mcf_pass_ncpi0, &b_mcf_pass_ncpi0);
   fChain->SetBranchAddress("mcf_pass_ccnopi", &mcf_pass_ccnopi, &b_mcf_pass_ccnopi);
   fChain->SetBranchAddress("mcf_pass_ncnopi", &mcf_pass_ncnopi, &b_mcf_pass_ncnopi);
   fChain->SetBranchAddress("mcf_pass_cccpi", &mcf_pass_cccpi, &b_mcf_pass_cccpi);
   fChain->SetBranchAddress("mcf_pass_nccpi", &mcf_pass_nccpi, &b_mcf_pass_nccpi);
   fChain->SetBranchAddress("X_SpcPts_v", &X_SpcPts_v, &b_X_SpcPts_v);
   fChain->SetBranchAddress("Y_SpcPts_v", &Y_SpcPts_v, &b_Y_SpcPts_v);
   fChain->SetBranchAddress("Z_SpcPts_v", &Z_SpcPts_v, &b_Z_SpcPts_v);
   fChain->SetBranchAddress("shr_id_MCStool", &shr_id_MCStool, &b_shr_pfp_id);
   fChain->SetBranchAddress("shr_hits_max_MCStool", &shr_hits_max_MCStool, &b_shr_hits_max_MCStool);
   fChain->SetBranchAddress("n_showers_contained_MCStool", &n_showers_contained_MCStool, &b_n_showers_contained_MCStool);
   fChain->SetBranchAddress("trkshrscore_v", &trkshrscore_v, &b_trkshrscore_v);
   fChain->SetBranchAddress("shrPCA_1Cr", &shrPCA_1Cr, &b_shrPCA_1Cr);
   fChain->SetBranchAddress("shrPCA_2Cr", &shrPCA_2Cr, &b_shrPCA_2Cr);
   fChain->SetBranchAddress("shrPCA_3Cr", &shrPCA_3Cr, &b_shrPCA_3Cr);
   fChain->SetBranchAddress("shrPCA_1Ce", &shrPCA_1Ce, &b_shrPCA_1Ce);
   fChain->SetBranchAddress("shrPCA_2Ce", &shrPCA_2Ce, &b_shrPCA_2Ce);
   fChain->SetBranchAddress("shrPCA_3Ce", &shrPCA_3Ce, &b_shrPCA_3Ce);
   fChain->SetBranchAddress("shrPCA1CAS", &shrPCA1CAS, &b_shrPCA1CAS);
   fChain->SetBranchAddress("shrPCA2CAS", &shrPCA2CAS, &b_shrPCA2CAS);
   fChain->SetBranchAddress("shrPCA3CAS", &shrPCA3CAS, &b_shrPCA3CAS);
   fChain->SetBranchAddress("shrPCA_1Cr2h", &shrPCA_1Cr2h, &b_shrPCA_1Cr2h);
   fChain->SetBranchAddress("shrPCA_2Cr2h", &shrPCA_2Cr2h, &b_shrPCA_2Cr2h);
   fChain->SetBranchAddress("shrPCA_3Cr2h", &shrPCA_3Cr2h, &b_shrPCA_3Cr2h);
   fChain->SetBranchAddress("shrPCA_1Cr1h", &shrPCA_1Cr1h, &b_shrPCA_1Cr1h);
   fChain->SetBranchAddress("shrPCA_2Cr1h", &shrPCA_2Cr1h, &b_shrPCA_2Cr1h);
   fChain->SetBranchAddress("shrPCA_3Cr1h", &shrPCA_3Cr1h, &b_shrPCA_3Cr1h);
   fChain->SetBranchAddress("shrMCSMom", &shrMCSMom, &b_shrMCSMom);
   fChain->SetBranchAddress("shrMCSMom1h", &shrMCSMom1h, &b_shrMCSMom1h);
   fChain->SetBranchAddress("shrMCSMom2h", &shrMCSMom2h, &b_shrMCSMom2h);
   fChain->SetBranchAddress("shrPCALen", &shrPCALen, &b_shrPCALen);
   fChain->SetBranchAddress("n_shrSpcPts", &n_shrSpcPts, &b_n_shrSpcPts);
   fChain->SetBranchAddress("PCAWin_1Cr_5cm", &PCAWin_1Cr_5cm, &b_PCAWin_1Cr_5cm);
   fChain->SetBranchAddress("PCAWin_2Cr_5cm", &PCAWin_2Cr_5cm, &b_PCAWin_2Cr_5cm);
   fChain->SetBranchAddress("PCAWin_3Cr_5cm", &PCAWin_3Cr_5cm, &b_PCAWin_3Cr_5cm);
   fChain->SetBranchAddress("PCAWin_dist_5cm", &PCAWin_dist_5cm, &b_PCAWin_dist_5cm);
   fChain->SetBranchAddress("PCAWin_npts_5cm", &PCAWin_npts_5cm, &b_PCAWin_npts_5cm);
   fChain->SetBranchAddress("shrStart_5cm", &shrStart_5cm, &b_shrStart_5cm);
   fChain->SetBranchAddress("shrStartMCS_5cm", &shrStartMCS_5cm, &b_shrStartMCS_5cm);
   fChain->SetBranchAddress("shrMCSAS_5cm", &shrMCSAS_5cm, &b_shrMCSAS_5cm);
   fChain->SetBranchAddress("shrPCA1CAS_5cm", &shrPCA1CAS_5cm, &b_shrPCA1CAS_5cm);
   fChain->SetBranchAddress("shrPCA2CAS_5cm", &shrPCA2CAS_5cm, &b_shrPCA2CAS_5cm);
   fChain->SetBranchAddress("shrPCA3CAS_5cm", &shrPCA3CAS_5cm, &b_shrPCA3CAS_5cm);
   fChain->SetBranchAddress("shrPCA1CMed_5cm", &shrPCA1CMed_5cm, &b__shrPCA1CMed_5cm);
   fChain->SetBranchAddress("PCAWin_1Cr_2_5cm", &PCAWin_1Cr_2_5cm, &b_PCAWin_1Cr_2_5cm);
   fChain->SetBranchAddress("PCAWin_2Cr_2_5cm", &PCAWin_2Cr_2_5cm, &b_PCAWin_2Cr_2_5cm);
   fChain->SetBranchAddress("PCAWin_3Cr_2_5cm", &PCAWin_3Cr_2_5cm, &b_PCAWin_3Cr_2_5cm);
   fChain->SetBranchAddress("PCAWin_dist_2_5cm", &PCAWin_dist_2_5cm, &b_PCAWin_dist_2_5cm);
   fChain->SetBranchAddress("PCAWin_npts_2_5cm", &PCAWin_npts_2_5cm, &b_PCAWin_npts_2_5cm);
   fChain->SetBranchAddress("shrStart_2_5cm", &shrStart_2_5cm, &b_shrStart_2_5cm);
   fChain->SetBranchAddress("shrStartMCS_2_5cm", &shrStartMCS_2_5cm, &b_shrStartMCS_2_5cm);
   fChain->SetBranchAddress("shrMCSAS_2_5cm", &shrMCSAS_2_5cm, &b_shrMCSAS_2_5cm);
   fChain->SetBranchAddress("shrPCA1CAS_2_5cm", &shrPCA1CAS_2_5cm, &b_shrPCA1CAS_2_5cm);
   fChain->SetBranchAddress("shrPCA2CAS_2_5cm", &shrPCA2CAS_2_5cm, &b_shrPCA2CAS_2_5cm);
   fChain->SetBranchAddress("shrPCA3CAS_2_5cm", &shrPCA3CAS_2_5cm, &b_shrPCA3CAS_2_5cm);
   fChain->SetBranchAddress("shrPCA1CMed_2_5cm", &shrPCA1CMed_2_5cm, &b__shrPCA1CMed_2_5cm);
   fChain->SetBranchAddress("DeltaMed", &DeltaMed, &b_DeltaMed);
   fChain->SetBranchAddress("DeltaMed1h", &DeltaMed1h, &b_DeltaMed1h);
   fChain->SetBranchAddress("DeltaMed2h", &DeltaMed2h, &b_DeltaMed2h);
   fChain->SetBranchAddress("DeltaRMS", &DeltaRMS, &b_DeltaRMS);
   fChain->SetBranchAddress("DeltaRMS1h", &DeltaRMS1h, &b_DeltaRMS1h);
   fChain->SetBranchAddress("DeltaRMS2h", &DeltaRMS2h, &b_DeltaRMS2h);
   fChain->SetBranchAddress("CylFrac_1cm", &CylFrac_1cm, &b_CylFrac_1cm);
   fChain->SetBranchAddress("CylFrac1h_1cm", &CylFrac1h_1cm, &b_CylFrac1h_1cm);
   fChain->SetBranchAddress("CylFrac2h_1cm", &CylFrac2h_1cm, &b_CylFrac2h_1cm);
   fChain->SetBranchAddress("CylFrac_2cm", &CylFrac_2cm, &b_CylFrac_2cm);
   fChain->SetBranchAddress("CylFrac1h_2cm", &CylFrac1h_2cm, &b_CylFrac1h_2cm);
   fChain->SetBranchAddress("CylFrac2h_2cm", &CylFrac2h_2cm, &b_CylFrac2h_2cm);
   fChain->SetBranchAddress("CylFrac_3cm", &CylFrac_3cm, &b_CylFrac_3cm);
   fChain->SetBranchAddress("CylFrac1h_3cm", &CylFrac1h_3cm, &b_CylFrac1h_3cm);
   fChain->SetBranchAddress("CylFrac2h_3cm", &CylFrac2h_3cm, &b_CylFrac2h_3cm);
   fChain->SetBranchAddress("CylFrac_4cm", &CylFrac_4cm, &b_CylFrac_4cm);
   fChain->SetBranchAddress("CylFrac1h_4cm", &CylFrac1h_4cm, &b_CylFrac1h_4cm);
   fChain->SetBranchAddress("CylFrac2h_4cm", &CylFrac2h_4cm, &b_CylFrac2h_4cm);
   fChain->SetBranchAddress("CylFrac_5cm", &CylFrac_5cm, &b_CylFrac_5cm);
   fChain->SetBranchAddress("CylFrac1h_5cm", &CylFrac1h_5cm, &b_CylFrac1h_5cm);
   fChain->SetBranchAddress("CylFrac2h_5cm", &CylFrac2h_5cm, &b_CylFrac2h_5cm);
   fChain->SetBranchAddress("NeutrinoEnergy0", &NeutrinoEnergy0, &b_NeutrinoEnergy0);
   fChain->SetBranchAddress("NeutrinoEnergy1", &NeutrinoEnergy1, &b_NeutrinoEnergy1);
   fChain->SetBranchAddress("NeutrinoEnergy2", &NeutrinoEnergy2, &b_NeutrinoEnergy2);
   fChain->SetBranchAddress("SliceCaloEnergy0", &SliceCaloEnergy0, &b_SliceCaloEnergy0);
   fChain->SetBranchAddress("SliceCaloEnergy1", &SliceCaloEnergy1, &b_SliceCaloEnergy1);
   fChain->SetBranchAddress("SliceCaloEnergy2", &SliceCaloEnergy2, &b_SliceCaloEnergy2);
   fChain->SetBranchAddress("pi0_mcgamma0_e", &pi0_mcgamma0_e, &b_pi0_mcgamma0_e);
   fChain->SetBranchAddress("pi0_mcgamma0_px", &pi0_mcgamma0_px, &b_pi0_mcgamma0_px);
   fChain->SetBranchAddress("pi0_mcgamma0_py", &pi0_mcgamma0_py, &b_pi0_mcgamma0_py);
   fChain->SetBranchAddress("pi0_mcgamma0_pz", &pi0_mcgamma0_pz, &b_pi0_mcgamma0_pz);
   fChain->SetBranchAddress("pi0_mcrcdot0", &pi0_mcrcdot0, &b_pi0_mcrcdot0);
   fChain->SetBranchAddress("pi0_mcrce0", &pi0_mcrce0, &b_pi0_mcrce0);
   fChain->SetBranchAddress("pi0_mcgamma1_e", &pi0_mcgamma1_e, &b_pi0_mcgamma1_e);
   fChain->SetBranchAddress("pi0_mcgamma1_px", &pi0_mcgamma1_px, &b_pi0_mcgamma1_px);
   fChain->SetBranchAddress("pi0_mcgamma1_py", &pi0_mcgamma1_py, &b_pi0_mcgamma1_py);
   fChain->SetBranchAddress("pi0_mcgamma1_pz", &pi0_mcgamma1_pz, &b_pi0_mcgamma1_pz);
   fChain->SetBranchAddress("pi0_mcrcdot1", &pi0_mcrcdot1, &b_pi0_mcrcdot1);
   fChain->SetBranchAddress("pi0_mcrce1", &pi0_mcrce1, &b_pi0_mcrce1);
   fChain->SetBranchAddress("pi0_nshower", &pi0_nshower, &b_pi0_nshower);
   fChain->SetBranchAddress("pi0_ntrack", &pi0_ntrack, &b_pi0_ntrack);
   fChain->SetBranchAddress("pi0_ngamma", &pi0_ngamma, &b_pi0_ngamma);
   fChain->SetBranchAddress("pi0_radlen1", &pi0_radlen1, &b_pi0_radlen1);
   fChain->SetBranchAddress("pi0_radlen2", &pi0_radlen2, &b_pi0_radlen2);
   fChain->SetBranchAddress("pi0_dot1", &pi0_dot1, &b_pi0_dot1);
   fChain->SetBranchAddress("pi0_dot2", &pi0_dot2, &b_pi0_dot2);
   fChain->SetBranchAddress("pi0_energy1_Y", &pi0_energy1_Y, &b_pi0_energy1_Y);
   fChain->SetBranchAddress("pi0_energy2_Y", &pi0_energy2_Y, &b_pi0_energy2_Y);
   fChain->SetBranchAddress("pi0_dir1_x", &pi0_dir1_x, &b_pi0_dir1_x);
   fChain->SetBranchAddress("pi0_dir1_y", &pi0_dir1_y, &b_pi0_dir1_y);
   fChain->SetBranchAddress("pi0_dir1_z", &pi0_dir1_z, &b_pi0_dir1_z);
   fChain->SetBranchAddress("pi0_dir2_x", &pi0_dir2_x, &b_pi0_dir2_x);
   fChain->SetBranchAddress("pi0_dir2_y", &pi0_dir2_y, &b_pi0_dir2_y);
   fChain->SetBranchAddress("pi0_dir2_z", &pi0_dir2_z, &b_pi0_dir2_z);
   fChain->SetBranchAddress("pi0_dedx1_Y", &pi0_dedx1_Y, &b_pi0_dedx1_Y);
   fChain->SetBranchAddress("pi0_dedx2_Y", &pi0_dedx2_Y, &b_pi0_dedx2_Y);
   fChain->SetBranchAddress("pi0_dedx1_fit_Y", &pi0_dedx1_fit_Y, &b_pi0_dedx1_fit_Y);
   fChain->SetBranchAddress("pi0_dedx2_fit_Y", &pi0_dedx2_fit_Y, &b_pi0_dedx2_fit_Y);
   fChain->SetBranchAddress("pi0_energy1_V", &pi0_energy1_V, &b_pi0_energy1_V);
   fChain->SetBranchAddress("pi0_energy2_V", &pi0_energy2_V, &b_pi0_energy2_V);
   fChain->SetBranchAddress("pi0_dedx1_V", &pi0_dedx1_V, &b_pi0_dedx1_V);
   fChain->SetBranchAddress("pi0_dedx2_V", &pi0_dedx2_V, &b_pi0_dedx2_V);
   fChain->SetBranchAddress("pi0_dedx1_fit_V", &pi0_dedx1_fit_V, &b_pi0_dedx1_fit_V);
   fChain->SetBranchAddress("pi0_dedx2_fit_V", &pi0_dedx2_fit_V, &b_pi0_dedx2_fit_V);
   fChain->SetBranchAddress("pi0_energy1_U", &pi0_energy1_U, &b_pi0_energy1_U);
   fChain->SetBranchAddress("pi0_energy2_U", &pi0_energy2_U, &b_pi0_energy2_U);
   fChain->SetBranchAddress("pi0_dedx1_U", &pi0_dedx1_U, &b_pi0_dedx1_U);
   fChain->SetBranchAddress("pi0_dedx2_U", &pi0_dedx2_U, &b_pi0_dedx2_U);
   fChain->SetBranchAddress("pi0_dedx1_fit_U", &pi0_dedx1_fit_U, &b_pi0_dedx1_fit_U);
   fChain->SetBranchAddress("pi0_dedx2_fit_U", &pi0_dedx2_fit_U, &b_pi0_dedx2_fit_U);
   fChain->SetBranchAddress("pi0_shrscore1", &pi0_shrscore1, &b_pi0_shrscore1);
   fChain->SetBranchAddress("pi0_shrscore2", &pi0_shrscore2, &b_pi0_shrscore2);
   fChain->SetBranchAddress("pi0_gammadot", &pi0_gammadot, &b_pi0_gammadot);
   fChain->SetBranchAddress("pi0_mass_Y", &pi0_mass_Y, &b_pi0_mass_Y);
   fChain->SetBranchAddress("pi0_mass_V", &pi0_mass_V, &b_pi0_mass_V);
   fChain->SetBranchAddress("pi0_mass_U", &pi0_mass_U, &b_pi0_mass_U);
   fChain->SetBranchAddress("pi0_rc_vtx_x", &pi0_rc_vtx_x, &b_pi0_rc_vtx_x);
   fChain->SetBranchAddress("pi0_rc_vtx_y", &pi0_rc_vtx_y, &b_pi0_rc_vtx_y);
   fChain->SetBranchAddress("pi0_rc_vtx_z", &pi0_rc_vtx_z, &b_pi0_rc_vtx_z);
   fChain->SetBranchAddress("pi0truth_gamma_parent", &pi0truth_gamma_parent, &b_pi0truth_gamma_parent);
   fChain->SetBranchAddress("pi0truth_elec_edep", &pi0truth_elec_edep, &b_pi0truth_elec_edep);
   fChain->SetBranchAddress("pi0truth_elec_etot", &pi0truth_elec_etot, &b_pi0truth_elec_etot);
   fChain->SetBranchAddress("pi0truth_elec_dist", &pi0truth_elec_dist, &b_pi0truth_elec_dist);
   fChain->SetBranchAddress("pi0truth_elec_parent", &pi0truth_elec_parent, &b_pi0truth_elec_parent);
   fChain->SetBranchAddress("pi0truth_gamma1_tid", &pi0truth_gamma1_tid, &b_pi0truth_gamma1_tid);
   fChain->SetBranchAddress("pi0truth_gamma1_edep", &pi0truth_gamma1_edep, &b_pi0truth_gamma1_edep);
   fChain->SetBranchAddress("pi0truth_gamma1_etot", &pi0truth_gamma1_etot, &b_pi0truth_gamma1_etot);
   fChain->SetBranchAddress("pi0truth_gamma1_dist", &pi0truth_gamma1_dist, &b_pi0truth_gamma1_dist);
   fChain->SetBranchAddress("pi0truth_gamma1_elec1", &pi0truth_gamma1_elec1, &b_pi0truth_gamma1_elec1);
   fChain->SetBranchAddress("pi0truth_gamma1_elec2", &pi0truth_gamma1_elec2, &b_pi0truth_gamma1_elec2);
   fChain->SetBranchAddress("pi0truth_gamma1_xpos", &pi0truth_gamma1_xpos, &b_pi0truth_gamma1_xpos);
   fChain->SetBranchAddress("pi0truth_gamma1_ypos", &pi0truth_gamma1_ypos, &b_pi0truth_gamma1_ypos);
   fChain->SetBranchAddress("pi0truth_gamma1_zpos", &pi0truth_gamma1_zpos, &b_pi0truth_gamma1_zpos);
   fChain->SetBranchAddress("pi0truth_gamma2_tid", &pi0truth_gamma2_tid, &b_pi0truth_gamma2_tid);
   fChain->SetBranchAddress("pi0truth_gamma2_edep", &pi0truth_gamma2_edep, &b_pi0truth_gamma2_edep);
   fChain->SetBranchAddress("pi0truth_gamma2_etot", &pi0truth_gamma2_etot, &b_pi0truth_gamma2_etot);
   fChain->SetBranchAddress("pi0truth_gamma2_dist", &pi0truth_gamma2_dist, &b_pi0truth_gamma2_dist);
   fChain->SetBranchAddress("pi0truth_gamma2_elec1", &pi0truth_gamma2_elec1, &b_pi0truth_gamma2_elec1);
   fChain->SetBranchAddress("pi0truth_gamma2_elec2", &pi0truth_gamma2_elec2, &b_pi0truth_gamma2_elec2);
   fChain->SetBranchAddress("pi0truth_gamma2_xpos", &pi0truth_gamma2_xpos, &b_pi0truth_gamma2_xpos);
   fChain->SetBranchAddress("pi0truth_gamma2_ypos", &pi0truth_gamma2_ypos, &b_pi0truth_gamma2_ypos);
   fChain->SetBranchAddress("pi0truth_gamma2_zpos", &pi0truth_gamma2_zpos, &b_pi0truth_gamma2_zpos);
   fChain->SetBranchAddress("pi0truth_gammadot", &pi0truth_gammadot, &b_pi0truth_gammadot);
   fChain->SetBranchAddress("pi0truth_run", &pi0truth_run, &b_pi0truth_run);
   fChain->SetBranchAddress("pi0truth_sub", &pi0truth_sub, &b_pi0truth_sub);
   fChain->SetBranchAddress("pi0truth_evt", &pi0truth_evt, &b_pi0truth_evt);
   fChain->SetBranchAddress("nflag_pl1", &nflag_pl1, &b_nflag_pl1);
   fChain->SetBranchAddress("nnoise_pl1", &nnoise_pl1, &b_nnoise_pl1);
   fChain->SetBranchAddress("nslhits_pl1", &nslhits_pl1, &b_nslhits_pl1);
   fChain->SetBranchAddress("nslnoise_pl1", &nslnoise_pl1, &b_nslnoise_pl1);
   fChain->SetBranchAddress("nhits_pl1", &nhits_pl1, &b_nhits_pl1);
   fChain->SetBranchAddress("frac_slnoise_pl1", &frac_slnoise_pl1, &b_frac_slnoise_pl1);
   fChain->SetBranchAddress("secondshower_U_charge", &secondshower_U_charge, &b_secondshower_U_charge);
   fChain->SetBranchAddress("secondshower_U_nhit", &secondshower_U_nhit, &b_secondshower_U_nhit);
   fChain->SetBranchAddress("secondshower_U_vtxdist", &secondshower_U_vtxdist, &b_secondshower_U_vtxdist);
   fChain->SetBranchAddress("secondshower_U_eigenratio", &secondshower_U_eigenratio, &b_secondshower_U_eigenratio);
   fChain->SetBranchAddress("secondshower_U_dot", &secondshower_U_dot, &b_secondshower_U_dot);
   fChain->SetBranchAddress("secondshower_U_dir", &secondshower_U_dir, &b_secondshower_U_dir);
   fChain->SetBranchAddress("secondshower_V_charge", &secondshower_V_charge, &b_secondshower_V_charge);
   fChain->SetBranchAddress("secondshower_V_nhit", &secondshower_V_nhit, &b_secondshower_V_nhit);
   fChain->SetBranchAddress("secondshower_V_vtxdist", &secondshower_V_vtxdist, &b_secondshower_V_vtxdist);
   fChain->SetBranchAddress("secondshower_V_eigenratio", &secondshower_V_eigenratio, &b_secondshower_V_eigenratio);
   fChain->SetBranchAddress("secondshower_V_dot", &secondshower_V_dot, &b_secondshower_V_dot);
   fChain->SetBranchAddress("secondshower_V_dir", &secondshower_V_dir, &b_secondshower_V_dir);
   fChain->SetBranchAddress("secondshower_Y_charge", &secondshower_Y_charge, &b_secondshower_Y_charge);
   fChain->SetBranchAddress("secondshower_Y_nhit", &secondshower_Y_nhit, &b_secondshower_Y_nhit);
   fChain->SetBranchAddress("secondshower_Y_vtxdist", &secondshower_Y_vtxdist, &b_secondshower_Y_vtxdist);
   fChain->SetBranchAddress("secondshower_Y_eigenratio", &secondshower_Y_eigenratio, &b_secondshower_Y_eigenratio);
   fChain->SetBranchAddress("secondshower_Y_dot", &secondshower_Y_dot, &b_secondshower_Y_dot);
   fChain->SetBranchAddress("secondshower_Y_dir", &secondshower_Y_dir, &b_secondshower_Y_dir);
   fChain->SetBranchAddress("shr_dedx_u_v", &shr_dedx_u_v, &b_shr_dedx_u_v);
   fChain->SetBranchAddress("shr_dedx_v_v", &shr_dedx_v_v, &b_shr_dedx_v_v);
   fChain->SetBranchAddress("shr_dedx_y_v", &shr_dedx_y_v, &b_shr_dedx_y_v);
   fChain->SetBranchAddress("shr_energy_u_v", &shr_energy_u_v, &b_shr_energy_u_v);
   fChain->SetBranchAddress("shr_energy_v_v", &shr_energy_v_v, &b_shr_energy_v_v);
   fChain->SetBranchAddress("shr_energy_y_v", &shr_energy_y_v, &b_shr_energy_y_v);
   fChain->SetBranchAddress("shr_pfp_id_v", &shr_pfp_id_v, &b_shr_pfp_id_v);
   fChain->SetBranchAddress("shr_start_x_v", &shr_start_x_v, &b_shr_start_x_v);
   fChain->SetBranchAddress("shr_start_y_v", &shr_start_y_v, &b_shr_start_y_v);
   fChain->SetBranchAddress("shr_start_z_v", &shr_start_z_v, &b_shr_start_z_v);
   fChain->SetBranchAddress("shr_dist_v", &shr_dist_v, &b_shr_dist_v);
   fChain->SetBranchAddress("shr_start_U_v", &shr_start_U_v, &b_shr_start_U_v);
   fChain->SetBranchAddress("shr_start_V_v", &shr_start_V_v, &b_shr_start_V_v);
   fChain->SetBranchAddress("shr_px_v", &shr_px_v, &b_shr_px_v);
   fChain->SetBranchAddress("shr_py_v", &shr_py_v, &b_shr_py_v);
   fChain->SetBranchAddress("shr_pz_v", &shr_pz_v, &b_shr_pz_v);
   fChain->SetBranchAddress("shr_openangle_v", &shr_openangle_v, &b_shr_openangle_v);
   fChain->SetBranchAddress("shr_theta_v", &shr_theta_v, &b_shr_theta_v);
   fChain->SetBranchAddress("shr_phi_v", &shr_phi_v, &b_shr_phi_v);
   fChain->SetBranchAddress("shr_pitch_u_v", &shr_pitch_u_v, &b_shr_pitch_u_v);
   fChain->SetBranchAddress("shr_pitch_v_v", &shr_pitch_v_v, &b_shr_pitch_v_v);
   fChain->SetBranchAddress("shr_pitch_y_v", &shr_pitch_y_v, &b_shr_pitch_y_v);
   fChain->SetBranchAddress("shr_tkfit_nhits_v", &shr_tkfit_nhits_v, &b_shr_tkfit_nhits_v);
   fChain->SetBranchAddress("shr_tkfit_start_x_v", &shr_tkfit_start_x_v, &b_shr_tkfit_start_x_v);
   fChain->SetBranchAddress("shr_tkfit_start_y_v", &shr_tkfit_start_y_v, &b_shr_tkfit_start_y_v);
   fChain->SetBranchAddress("shr_tkfit_start_z_v", &shr_tkfit_start_z_v, &b_shr_tkfit_start_z_v);
   fChain->SetBranchAddress("shr_tkfit_start_U_v", &shr_tkfit_start_U_v, &b_shr_tkfit_start_U_v);
   fChain->SetBranchAddress("shr_tkfit_start_V_v", &shr_tkfit_start_V_v, &b_shr_tkfit_start_V_v);
   fChain->SetBranchAddress("shr_tkfit_theta_v", &shr_tkfit_theta_v, &b_shr_tkfit_theta_v);
   fChain->SetBranchAddress("shr_tkfit_phi_v", &shr_tkfit_phi_v, &b_shr_tkfit_phi_v);
   fChain->SetBranchAddress("shr_tkfit_pitch_u_v", &shr_tkfit_pitch_u_v, &b_shr_tkfit_pitch_u_v);
   fChain->SetBranchAddress("shr_tkfit_pitch_v_v", &shr_tkfit_pitch_v_v, &b_shr_tkfit_pitch_v_v);
   fChain->SetBranchAddress("shr_tkfit_pitch_y_v", &shr_tkfit_pitch_y_v, &b_shr_tkfit_pitch_y_v);
   fChain->SetBranchAddress("shr_tkfit_dedx_u_v", &shr_tkfit_dedx_u_v, &b_shr_tkfit_dedx_u_v);
   fChain->SetBranchAddress("shr_tkfit_dedx_v_v", &shr_tkfit_dedx_v_v, &b_shr_tkfit_dedx_v_v);
   fChain->SetBranchAddress("shr_tkfit_dedx_y_v", &shr_tkfit_dedx_y_v, &b_shr_tkfit_dedx_y_v);
   fChain->SetBranchAddress("shr_tkfit_gap10_dedx_u_v", &shr_tkfit_gap10_dedx_u_v, &b_shr_tkfit_gap10_dedx_u_v);
   fChain->SetBranchAddress("shr_tkfit_gap10_dedx_v_v", &shr_tkfit_gap10_dedx_v_v, &b_shr_tkfit_gap10_dedx_v_v);
   fChain->SetBranchAddress("shr_tkfit_gap10_dedx_y_v", &shr_tkfit_gap10_dedx_y_v, &b_shr_tkfit_gap10_dedx_y_v);
   fChain->SetBranchAddress("shr_tkfit_dedx_nhits_u_v", &shr_tkfit_dedx_nhits_u_v, &b_shr_tkfit_dedx_nhits_u_v);
   fChain->SetBranchAddress("shr_tkfit_dedx_nhits_v_v", &shr_tkfit_dedx_nhits_v_v, &b_shr_tkfit_dedx_nhits_v_v);
   fChain->SetBranchAddress("shr_tkfit_dedx_nhits_y_v", &shr_tkfit_dedx_nhits_y_v, &b_shr_tkfit_dedx_nhits_y_v);
   fChain->SetBranchAddress("shr_llr_pid_u_v", &shr_llr_pid_u_v, &b_shr_llr_pid_u_v);
   fChain->SetBranchAddress("shr_llr_pid_v_v", &shr_llr_pid_v_v, &b_shr_llr_pid_v_v);
   fChain->SetBranchAddress("shr_llr_pid_y_v", &shr_llr_pid_y_v, &b_shr_llr_pid_y_v);
   fChain->SetBranchAddress("shr_llr_pid_v", &shr_llr_pid_v, &b_shr_llr_pid_v);
   fChain->SetBranchAddress("shr_llr_pid_score_v", &shr_llr_pid_score_v, &b_shr_llr_pid_score_v);
   fChain->SetBranchAddress("shr_moliere_avg_v", &shr_moliere_avg_v, &b_shr_moliere_avg_v);
   fChain->SetBranchAddress("shr_moliere_rms_v", &shr_moliere_rms_v, &b_shr_moliere_rms_v);
   fChain->SetBranchAddress("origevnunhits", &origevnunhits, &b_origevnunhits);
   fChain->SetBranchAddress("evnunhits", &evnunhits, &b_evnunhits);
   fChain->SetBranchAddress("evlepnhits", &evlepnhits, &b_evlepnhits);
   fChain->SetBranchAddress("evpronhits", &evpronhits, &b_evpronhits);
   fChain->SetBranchAddress("evpi1nhits", &evpi1nhits, &b_evpi1nhits);
   fChain->SetBranchAddress("evpi0nhits", &evpi0nhits, &b_evpi0nhits);
   fChain->SetBranchAddress("evneunhits", &evneunhits, &b_evneunhits);
   fChain->SetBranchAddress("evgamnhits", &evgamnhits, &b_evgamnhits);
   fChain->SetBranchAddress("evothnhits", &evothnhits, &b_evothnhits);
   fChain->SetBranchAddress("slnunhits", &slnunhits, &b_slnunhits);
   fChain->SetBranchAddress("sllepnhits", &sllepnhits, &b_sllepnhits);
   fChain->SetBranchAddress("slpronhits", &slpronhits, &b_slpronhits);
   fChain->SetBranchAddress("slpi1nhits", &slpi1nhits, &b_slpi1nhits);
   fChain->SetBranchAddress("slpi0nhits", &slpi0nhits, &b_slpi0nhits);
   fChain->SetBranchAddress("slneunhits", &slneunhits, &b_slneunhits);
   fChain->SetBranchAddress("slgamnhits", &slgamnhits, &b_slgamnhits);
   fChain->SetBranchAddress("slothnhits", &slothnhits, &b_slothnhits);
   fChain->SetBranchAddress("pfnunhits", &pfnunhits, &b_pfnunhits);
   fChain->SetBranchAddress("pflepnhits", &pflepnhits, &b_pflepnhits);
   fChain->SetBranchAddress("pfpronhits", &pfpronhits, &b_pfpronhits);
   fChain->SetBranchAddress("pfpi1nhits", &pfpi1nhits, &b_pfpi1nhits);
   fChain->SetBranchAddress("pfpi0nhits", &pfpi0nhits, &b_pfpi0nhits);
   fChain->SetBranchAddress("pfneunhits", &pfneunhits, &b_pfneunhits);
   fChain->SetBranchAddress("pfgamnhits", &pfgamnhits, &b_pfgamnhits);
   fChain->SetBranchAddress("pfothnhits", &pfothnhits, &b_pfothnhits);
   fChain->SetBranchAddress("nu_completeness_from_pfp", &nu_completeness_from_pfp, &b_nu_completeness_from_pfp);
   fChain->SetBranchAddress("nu_purity_from_pfp", &nu_purity_from_pfp, &b_nu_purity_from_pfp);
   fChain->SetBranchAddress("trk_bragg_p_v", &trk_bragg_p_v, &b_trk_bragg_p_v);
   fChain->SetBranchAddress("trk_bragg_mu_v", &trk_bragg_mu_v, &b_trk_bragg_mu_v);
   fChain->SetBranchAddress("trk_bragg_pion_v", &trk_bragg_pion_v, &b_trk_bragg_pion_v);
   fChain->SetBranchAddress("trk_bragg_mip_v", &trk_bragg_mip_v, &b_trk_bragg_mip_v);
   fChain->SetBranchAddress("trk_bragg_p_alt_dir_v", &trk_bragg_p_alt_dir_v, &b_trk_bragg_p_alt_dir_v);
   fChain->SetBranchAddress("trk_bragg_mu_alt_dir_v", &trk_bragg_mu_alt_dir_v, &b_trk_bragg_mu_alt_dir_v);
   fChain->SetBranchAddress("trk_bragg_pion_alt_dir_v", &trk_bragg_pion_alt_dir_v, &b_trk_bragg_pion_alt_dir_v);
   fChain->SetBranchAddress("trk_bragg_p_fwd_preferred_v", &trk_bragg_p_fwd_preferred_v, &b_trk_bragg_p_fwd_preferred_v);
   fChain->SetBranchAddress("trk_bragg_mu_fwd_preferred_v", &trk_bragg_mu_fwd_preferred_v, &b_trk_bragg_mu_fwd_preferred_v);
   fChain->SetBranchAddress("trk_bragg_pion_fwd_preferred_v", &trk_bragg_pion_fwd_preferred_v, &b_trk_bragg_pion_fwd_preferred_v);
   fChain->SetBranchAddress("trk_pida_v", &trk_pida_v, &b_trk_pida_v);
   fChain->SetBranchAddress("trk_pid_chipr_v", &trk_pid_chipr_v, &b_trk_pid_chipr_v);
   fChain->SetBranchAddress("trk_pid_chipi_v", &trk_pid_chipi_v, &b_trk_pid_chipi_v);
   fChain->SetBranchAddress("trk_pid_chika_v", &trk_pid_chika_v, &b_trk_pid_chika_v);
   fChain->SetBranchAddress("trk_pid_chimu_v", &trk_pid_chimu_v, &b_trk_pid_chimu_v);
   fChain->SetBranchAddress("trk_bragg_p_u_v", &trk_bragg_p_u_v, &b_trk_bragg_p_u_v);
   fChain->SetBranchAddress("trk_bragg_mu_u_v", &trk_bragg_mu_u_v, &b_trk_bragg_mu_u_v);
   fChain->SetBranchAddress("trk_bragg_pion_u_v", &trk_bragg_pion_u_v, &b_trk_bragg_pion_u_v);
   fChain->SetBranchAddress("trk_bragg_mip_u_v", &trk_bragg_mip_u_v, &b_trk_bragg_mip_u_v);
   fChain->SetBranchAddress("trk_bragg_p_alt_dir_u_v", &trk_bragg_p_alt_dir_u_v, &b_trk_bragg_p_alt_dir_u_v);
   fChain->SetBranchAddress("trk_bragg_mu_alt_dir_u_v", &trk_bragg_mu_alt_dir_u_v, &b_trk_bragg_mu_alt_dir_u_v);
   fChain->SetBranchAddress("trk_bragg_pion_alt_dir_u_v", &trk_bragg_pion_alt_dir_u_v, &b_trk_bragg_pion_alt_dir_u_v);
   fChain->SetBranchAddress("trk_bragg_p_fwd_preferred_u_v", &trk_bragg_p_fwd_preferred_u_v, &b_trk_bragg_p_fwd_preferred_u_v);
   fChain->SetBranchAddress("trk_bragg_mu_fwd_preferred_u_v", &trk_bragg_mu_fwd_preferred_u_v, &b_trk_bragg_mu_fwd_preferred_u_v);
   fChain->SetBranchAddress("trk_bragg_pion_fwd_preferred_u_v", &trk_bragg_pion_fwd_preferred_u_v, &b_trk_bragg_pion_fwd_preferred_u_v);
   fChain->SetBranchAddress("trk_pida_u_v", &trk_pida_u_v, &b_trk_pida_u_v);
   fChain->SetBranchAddress("trk_pid_chipr_u_v", &trk_pid_chipr_u_v, &b_trk_pid_chipr_u_v);
   fChain->SetBranchAddress("trk_pid_chipi_u_v", &trk_pid_chipi_u_v, &b_trk_pid_chipi_u_v);
   fChain->SetBranchAddress("trk_pid_chika_u_v", &trk_pid_chika_u_v, &b_trk_pid_chika_u_v);
   fChain->SetBranchAddress("trk_pid_chimu_u_v", &trk_pid_chimu_u_v, &b_trk_pid_chimu_u_v);
   fChain->SetBranchAddress("trk_bragg_p_v_v", &trk_bragg_p_v_v, &b_trk_bragg_p_v_v);
   fChain->SetBranchAddress("trk_bragg_mu_v_v", &trk_bragg_mu_v_v, &b_trk_bragg_mu_v_v);
   fChain->SetBranchAddress("trk_bragg_pion_v_v", &trk_bragg_pion_v_v, &b_trk_bragg_pion_v_v);
   fChain->SetBranchAddress("trk_bragg_mip_v_v", &trk_bragg_mip_v_v, &b_trk_bragg_mip_v_v);
   fChain->SetBranchAddress("trk_bragg_p_alt_dir_v_v", &trk_bragg_p_alt_dir_v_v, &b_trk_bragg_p_alt_dir_v_v);
   fChain->SetBranchAddress("trk_bragg_mu_alt_dir_v_v", &trk_bragg_mu_alt_dir_v_v, &b_trk_bragg_mu_alt_dir_v_v);
   fChain->SetBranchAddress("trk_bragg_pion_alt_dir_v_v", &trk_bragg_pion_alt_dir_v_v, &b_trk_bragg_pion_alt_dir_v_v);
   fChain->SetBranchAddress("trk_bragg_p_fwd_preferred_v_v", &trk_bragg_p_fwd_preferred_v_v, &b_trk_bragg_p_fwd_preferred_v_v);
   fChain->SetBranchAddress("trk_bragg_mu_fwd_preferred_v_v", &trk_bragg_mu_fwd_preferred_v_v, &b_trk_bragg_mu_fwd_preferred_v_v);
   fChain->SetBranchAddress("trk_bragg_pion_fwd_preferred_v_v", &trk_bragg_pion_fwd_preferred_v_v, &b_trk_bragg_pion_fwd_preferred_v_v);
   fChain->SetBranchAddress("trk_pida_v_v", &trk_pida_v_v, &b_trk_pida_v_v);
   fChain->SetBranchAddress("trk_pid_chipr_v_v", &trk_pid_chipr_v_v, &b_trk_pid_chipr_v_v);
   fChain->SetBranchAddress("trk_pid_chipi_v_v", &trk_pid_chipi_v_v, &b_trk_pid_chipi_v_v);
   fChain->SetBranchAddress("trk_pid_chika_v_v", &trk_pid_chika_v_v, &b_trk_pid_chika_v_v);
   fChain->SetBranchAddress("trk_pid_chimu_v_v", &trk_pid_chimu_v_v, &b_trk_pid_chimu_v_v);
   fChain->SetBranchAddress("trk_pfp_id_v", &trk_pfp_id_v, &b_trk_pfp_id_v);
   fChain->SetBranchAddress("trk_dir_x_v", &trk_dir_x_v, &b_trk_dir_x_v);
   fChain->SetBranchAddress("trk_dir_y_v", &trk_dir_y_v, &b_trk_dir_y_v);
   fChain->SetBranchAddress("trk_dir_z_v", &trk_dir_z_v, &b_trk_dir_z_v);
   fChain->SetBranchAddress("trk_start_x_v", &trk_start_x_v, &b_trk_start_x_v);
   fChain->SetBranchAddress("trk_start_y_v", &trk_start_y_v, &b_trk_start_y_v);
   fChain->SetBranchAddress("trk_start_z_v", &trk_start_z_v, &b_trk_start_z_v);
   fChain->SetBranchAddress("trk_sce_start_x_v", &trk_sce_start_x_v, &b_trk_sce_start_x_v);
   fChain->SetBranchAddress("trk_sce_start_y_v", &trk_sce_start_y_v, &b_trk_sce_start_y_v);
   fChain->SetBranchAddress("trk_sce_start_z_v", &trk_sce_start_z_v, &b_trk_sce_start_z_v);
   fChain->SetBranchAddress("trk_end_x_v", &trk_end_x_v, &b_trk_end_x_v);
   fChain->SetBranchAddress("trk_end_y_v", &trk_end_y_v, &b_trk_end_y_v);
   fChain->SetBranchAddress("trk_end_z_v", &trk_end_z_v, &b_trk_end_z_v);
   fChain->SetBranchAddress("trk_sce_end_x_v", &trk_sce_end_x_v, &b_trk_sce_end_x_v);
   fChain->SetBranchAddress("trk_sce_end_y_v", &trk_sce_end_y_v, &b_trk_sce_end_y_v);
   fChain->SetBranchAddress("trk_sce_end_z_v", &trk_sce_end_z_v, &b_trk_sce_end_z_v);
   fChain->SetBranchAddress("trk_distance_v", &trk_distance_v, &b_trk_distance_v);
   fChain->SetBranchAddress("trk_theta_v", &trk_theta_v, &b_trk_theta_v);
   fChain->SetBranchAddress("trk_phi_v", &trk_phi_v, &b_trk_phi_v);
   fChain->SetBranchAddress("trk_len_v", &trk_len_v, &b_trk_len_v);
   fChain->SetBranchAddress("trk_mcs_muon_mom_v", &trk_mcs_muon_mom_v, &b_trk_mcs_muon_mom_v);
   fChain->SetBranchAddress("trk_range_muon_mom_v", &trk_range_muon_mom_v, &b_trk_range_muon_mom_v);
   fChain->SetBranchAddress("trk_energy_proton_v", &trk_energy_proton_v, &b_trk_energy_proton_v);
   fChain->SetBranchAddress("trk_energy_muon_v", &trk_energy_muon_v, &b_trk_energy_muon_v);
   fChain->SetBranchAddress("trk_calo_energy_u_v", &trk_calo_energy_u_v, &b_trk_calo_energy_u_v);
   fChain->SetBranchAddress("trk_calo_energy_v_v", &trk_calo_energy_v_v, &b_trk_calo_energy_v_v);
   fChain->SetBranchAddress("trk_calo_energy_y_v", &trk_calo_energy_y_v, &b_trk_calo_energy_y_v);
   fChain->SetBranchAddress("trk_llr_pid_u_v", &trk_llr_pid_u_v, &b_trk_llr_pid_u_v);
   fChain->SetBranchAddress("trk_llr_pid_v_v", &trk_llr_pid_v_v, &b_trk_llr_pid_v_v);
   fChain->SetBranchAddress("trk_llr_pid_y_v", &trk_llr_pid_y_v, &b_trk_llr_pid_y_v);
   fChain->SetBranchAddress("trk_llr_pid_v", &trk_llr_pid_v, &b_trk_llr_pid_v);
   fChain->SetBranchAddress("trk_llr_pid_score_v", &trk_llr_pid_score_v, &b_trk_llr_pid_score_v);
   fChain->SetBranchAddress("trk_trunk_dEdx_u_v", &trk_trunk_dEdx_u_v, &b_trk_trunk_dEdx_u_v);
   fChain->SetBranchAddress("trk_trunk_dEdx_v_v", &trk_trunk_dEdx_v_v, &b_trk_trunk_dEdx_v_v);
   fChain->SetBranchAddress("trk_trunk_dEdx_y_v", &trk_trunk_dEdx_y_v, &b_trk_trunk_dEdx_y_v);
   fChain->SetBranchAddress("trk_trunk_rr_dEdx_u_v", &trk_trunk_rr_dEdx_u_v, &b_trk_trunk_rr_dEdx_u_v);
   fChain->SetBranchAddress("trk_trunk_rr_dEdx_v_v", &trk_trunk_rr_dEdx_v_v, &b_trk_trunk_rr_dEdx_v_v);
   fChain->SetBranchAddress("trk_trunk_rr_dEdx_y_v", &trk_trunk_rr_dEdx_y_v, &b_trk_trunk_rr_dEdx_y_v);
   fChain->SetBranchAddress("trk_nhits_u_v", &trk_nhits_u_v, &b_trk_nhits_u_v);
   fChain->SetBranchAddress("trk_nhits_v_v", &trk_nhits_v_v, &b_trk_nhits_v_v);
   fChain->SetBranchAddress("trk_nhits_y_v", &trk_nhits_y_v, &b_trk_nhits_y_v);
   fChain->SetBranchAddress("trk_avg_deflection_mean_v", &trk_avg_deflection_mean_v, &b_trk_avg_deflection_mean_v);
   fChain->SetBranchAddress("trk_avg_deflection_stdev_v", &trk_avg_deflection_stdev_v, &b_trk_avg_deflection_stdev_v);
   fChain->SetBranchAddress("trk_avg_deflection_separation_mean_v", &trk_avg_deflection_separation_mean_v, &b_trk_avg_deflection_separation_mean_v);
   fChain->SetBranchAddress("trk_end_spacepoints_v", &trk_end_spacepoints_v, &b_trk_end_spacepoints_v);
   fChain->SetBranchAddress("bdt_nuNCpi0", &bdt_nuNCpi0, &b_bdt_nuNCpi0);
   fChain->SetBranchAddress("bdt_numuCCpi0", &bdt_numuCCpi0, &b_bdt_numuCCpi0);
   fChain->SetBranchAddress("bdt_numuCC", &bdt_numuCC, &b_bdt_numuCC);
   fChain->SetBranchAddress("bdt_ext", &bdt_ext, &b_bdt_ext);
   fChain->SetBranchAddress("bdt_cosmic", &bdt_cosmic, &b_bdt_cosmic);
   fChain->SetBranchAddress("bdt_global", &bdt_global, &b_bdt_global);
   fChain->SetBranchAddress("pass_antibdt_filter", &pass_antibdt_filter, &b_bdt_global);
   fChain->SetBranchAddress("bdt_pi0_np", &bdt_pi0_np, &b_bdt_pi0_np);
   fChain->SetBranchAddress("bdt_nonpi0_np", &bdt_nonpi0_np, &b_bdt_nonpi0_np);
   fChain->SetBranchAddress("bdt_bkg_0p", &bdt_bkg_0p, &b_bdt_bkg_0p);
   fChain->SetBranchAddress("anglediff_Y", &anglediff_Y, &b_anglediff_Y);
   fChain->SetBranchAddress("anglediff_V", &anglediff_V, &b_anglediff_V);
   fChain->SetBranchAddress("anglediff_U", &anglediff_U, &b_anglediff_U);
   fChain->SetBranchAddress("trkpid", &trkpid, &b_trkpid);
   Notify();
}

bool NeutrinoSelectionFilter::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return true;
}

void NeutrinoSelectionFilter::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t NeutrinoSelectionFilter::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef NeutrinoSelectionFilter_cxx
