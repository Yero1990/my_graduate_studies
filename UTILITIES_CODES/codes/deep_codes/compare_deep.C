#include "../plotting_utilities/histUtils.C"   //load utilities code for histo plotting

//Script to make comparison between SIMC and Commissioning Data from HallC Spring 2018
//The comparisons are all charge normalized and corrected for all inefficiencies

void compare_deep(int pm, int set, int thnq)
{

  //gROOT->SetBatch(kTRUE);  
  //gStyle->SetOptStat(1001111);
  gStyle->SetOptStat(0000000);
  string thnq_set;
  TString root_dir;
  if(thnq==-1){thnq_set = "all_thnq_Q2_4to5GeV";}
  else{thnq_set = Form("thnq%d_Q2_4to5GeV", thnq);}

  
  if(thnq==-1){root_dir = Form("../../root_files/pm%d_fsiXsec_set%d_%s/", pm, set, thnq_set.c_str());}
  else{TString root_dir = Form("../../root_files/pm%d_fsiXsec_set%d_%s/", pm, set, thnq_set.c_str());}
  string simc_fname = Form("deep_simc_histos_pm%d_lagetfsi_rad_set%d.root", pm, set);
  string data_fname = Form("deep_data_histos_pm%d_set%d_combined.root", pm, set);
    
  TString simc_filename_fsi = root_dir+simc_fname;
  TString data_filename = root_dir+data_fname;
  cout << simc_filename_fsi << endl;
  cout << data_filename << endl;
  //Open SIMC/data ROOT files;
  TFile *simc_file_fsi = new TFile(simc_filename_fsi);
  TFile *data_file = new TFile(data_filename);


  //-------------Spectrometer Acceptance Plots--------
  
  //FSI
  TH1F *simc_eytar_fsi =  0;
  TH1F *simc_exptar_fsi =  0;
  TH1F *simc_eyptar_fsi =  0;
  TH1F *simc_edelta_fsi =  0;

  TH1F *simc_hytar_fsi =  0;
  TH1F *simc_hxptar_fsi =  0;
  TH1F *simc_hyptar_fsi =  0;
  TH1F *simc_hdelta_fsi =  0;

  
  //Define data histos
  TH1F *data_eytar = 0;
  TH1F *data_exptar =  0;
  TH1F *data_eyptar =  0;
  TH1F *data_edelta =  0;

  TH1F *data_hytar = 0;
  TH1F *data_hxptar =  0;
  TH1F *data_hyptar =  0;
  TH1F *data_hdelta =  0;


  //------------------ANALYSIS PLOTS------------------

  //SIMC Analysis Cuts (FSI)
  TH1F *simc_Q2_fsi =  0;
  TH1F *simc_thnq_fsi = 0;
  TH1F *simc_emiss_fsi = 0;
  TH1F *simc_ztar_diff_fsi = 0; 
  TH2F *simc_HMS_Coll_fsi = 0;


  //Data Analysis Plots
  TH1F *data_Q2 =  0;
  TH1F *data_thnq = 0;
  TH1F *data_emiss = 0;
  TH1F *data_ztar_diff = 0;
  TH2F *data_HMS_Coll = 0;
  
  TH1F *data_CoinTime = 0;
  TH1F *data_pid_eCal = 0;


  //Missing Momentum histograms 
  TH1F *data_Pm = 0;
  TH1F *simc_Pm_fsi = 0;
  

  //---------Get SIMC Histograms------

  //change to FSI simc_file
  simc_file_fsi->cd();
  
  simc_file_fsi->GetObject("H_ztar_diff_sys", simc_ztar_diff_fsi);
  simc_file_fsi->GetObject("H_hXColl_vs_hYColl_sys", simc_HMS_Coll_fsi);

  //------SIMC Spectrometer Acceptance------
  simc_file_fsi->GetObject("H_eytar", simc_eytar_fsi);
  simc_file_fsi->GetObject("H_exptar", simc_exptar_fsi);
  simc_file_fsi->GetObject("H_eyptar", simc_eyptar_fsi);
  simc_file_fsi->GetObject("H_edelta_sys", simc_edelta_fsi);

  simc_file_fsi->GetObject("H_hytar", simc_hytar_fsi);
  simc_file_fsi->GetObject("H_hxptar", simc_hxptar_fsi);
  simc_file_fsi->GetObject("H_hyptar", simc_hyptar_fsi);
  simc_file_fsi->GetObject("H_hdelta_sys", simc_hdelta_fsi);


  //------SIMC ANALYSIS PLOTS-----
  simc_file_fsi->GetObject("H_Q2_sys", simc_Q2_fsi);
  simc_file_fsi->GetObject("H_theta_nq_sys", simc_thnq_fsi);
  simc_file_fsi->GetObject("H_Em_nuc_sys", simc_emiss_fsi);
  simc_file_fsi->GetObject("H_Pm", simc_Pm_fsi);

 
  //Get Histogram objects from DATA rootfile
  //change to data_file
  data_file->cd();
  
  data_file->GetObject("H_ztar_diff_sys", data_ztar_diff);
  data_file->GetObject("H_hXColl_vs_hYColl_sys", data_HMS_Coll); 
  data_file->GetObject("H_pcal_etotTrkNorm_sys", data_pid_eCal); 
  data_file->GetObject("H_ctime_sys", data_CoinTime); 

  //-------DATA Spectrometer Acceptance-----
  data_file->GetObject("H_eytar", data_eytar);
  data_file->GetObject("H_exptar", data_exptar);
  data_file->GetObject("H_eyptar", data_eyptar);
  data_file->GetObject("H_edelta_sys", data_edelta);
  
  data_file->GetObject("H_hytar", data_hytar);
  data_file->GetObject("H_hxptar", data_hxptar);
  data_file->GetObject("H_hyptar", data_hyptar);
  data_file->GetObject("H_hdelta_sys", data_hdelta);

  //-----Get Data Analysis Histograms
  data_file->GetObject("H_Q2_sys", data_Q2);
  data_file->GetObject("H_theta_nq_sys", data_thnq);
  data_file->GetObject("H_Em_nuc_sys", data_emiss);  
  data_file->GetObject("H_Pm", data_Pm);

  
  //Plot the SHMS Recon
  //hist_ratio(data_eytar, simc_eytar_fsi, "SHMS Y_{tar} [cm]", "Counts", "SHMS Y_{tar}");
  //hist_ratio(data_exptar, simc_exptar_fsi, "SHMS X'_{tar} [rad]", "Counts", "SHMS X'_{tar}");
  //hist_ratio(data_eyptar, simc_eyptar_fsi, "SHMS Y'_{tar} [rad]", "Counts", "SHMS Y'_{tar}");
  //hist_ratio(data_edelta, simc_edelta_fsi, "SHMS #delta [%]", "Counts", "SHMS #delta");

  //Plot the HMS Recon
  //hist_ratio(data_hytar, simc_hytar_fsi, "HMS Y_{tar} [cm]", "Counts", "HMS Y_{tar}");
  //hist_ratio(data_hxptar, simc_hxptar_fsi, "HMS X'_{tar} [rad]", "Counts", "HMS X'_{tar}");
  //hist_ratio(data_hyptar, simc_hyptar_fsi, "HMS Y'_{tar} [rad]", "Counts", "HMS Y'_{tar}");
  //hist_ratio(data_hdelta, simc_hdelta_fsi, "HMS #delta [%]", "Counts", "HMS #delta");

  //Plot Only Data Histos
  //plot_hist(data_pid_eCal, "SHMS Calorimeter E_{dep}/P_{trk}", "Counts", "SHMS Calorimeter Total Normalized Energy", "logy");
  //plot_hist(data_CoinTime, "Coincidence Time [ns]", "Counts", "Coincidence Time", "logy");


  //ONLY PLOT DATA-2-SIMC COMPARISONS (NOT RATIOS)
  //compare_hist(data_hdelta, simc_hdelta_fsi, "HMS #delta [%]", "Counts", "HMS #delta");
  // compare_hist(data_edelta, simc_edelta_fsi, "SHMS #delta [%]", "Counts", "SHMS #delta");

  //compare_hist(data_Q2, simc_Q2_fsi, "Q^{2} [GeV^{2}]", "Counts", "4-Momentum Transfer, Q^{2}");
  //compare_hist(data_thnq, simc_thnq_fsi, "#theta_{nq} [deg]", "Counts", "Neutron Recoil Angles, #theta_{nq}");
  //compare_hist(data_ztar_diff, simc_ztar_diff_fsi, "Z_{tar} Difference [cm]", "Counts", "Z_{tar} Difference");
  //compare_hist(data_emiss, simc_emiss_fsi, "Missing Energy, E_{m} [GeV]", "Counts", "Nuclear Missing Energy");

  //compare_hist(data_Pm, simc_Pm_fsi, "Neutron Recoil Momentum, P_{r} [GeV]", "Counts", "Neutron Recoil Momentum");

    //ONLY COMPARE HISTOS  
   //Plot the SHMS Recon
  //compare_hist(data_eytar, simc_eytar_fsi, "SHMS Y_{tar} [cm]", "Counts", "SHMS Y_{tar}");
  //compare_hist(data_exptar, simc_exptar_fsi, "SHMS X'_{tar} [rad]", "Counts", "SHMS X'_{tar}");
  //compare_hist(data_eyptar, simc_eyptar_fsi, "SHMS Y'_{tar} [rad]", "Counts", "SHMS Y'_{tar}");
  //compare_hist(data_edelta, simc_edelta_fsi, "SHMS #delta [%]", "Counts", "SHMS #delta");

  //Plot the HMS Recon
  //compare_hist(data_hytar, simc_hytar_fsi, "HMS Y_{tar} [cm]", "Counts", "HMS Y_{tar}");
  //compare_hist(data_hxptar, simc_hxptar_fsi, "HMS X'_{tar} [rad]", "Counts", "HMS X'_{tar}");
  //compare_hist(data_hyptar, simc_hyptar_fsi, "HMS Y'_{tar} [rad]", "Counts", "HMS Y'_{tar}");
  //compare_hist(data_hdelta, simc_hdelta_fsi, "HMS #delta [%]", "Counts", "HMS #delta");
}
