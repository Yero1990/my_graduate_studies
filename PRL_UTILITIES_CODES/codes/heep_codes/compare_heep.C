#include "../plotting_utilities/histUtils.C"   //load utilities code for histo plotting

//Script to make comparison between SIMC and Commissioning Data from HallC Spring 2018
//The comparisons are all charge normalized and corrected for all inefficiencies

void compare_heep(int run)
{

  //gROOT->SetBatch(kTRUE);  
  gStyle->SetOptStat(0000000);


  TString root_dir = "../../root_files/pre_HEEP_ELASTICS/";
  //TString root_dir = "../../root_files/HEEP_ELASTICS/";

  //NO CUT and Ztar Diff Cut
  //TString data_filename = "../../root_files/heep_data_histos_3288_combined.root";    
  //TString simc_filename = "../../root_files/heep_simc_histos_3288_rad.root";       

  
  //Q2 4,5 GeV2 CUT
  //TString data_filename = "heep_data_histos_3288_combined_Q2_4to5.root";    
  //TString simc_filename = "heep_simc_histos_3288_rad_Q2_4to5.root";       

  //Full Q2 range
  //TString data_filename = Form("heep_data_histos_%d_combined_Q2full.root", run);    
  //TString simc_filename = Form("heep_simc_histos_%d_rad_Q2full.root", run);       

  
  //Pre-defined SIMC/data root file names containing histogram object to comapare
  TString simc_filename =  Form("heep_simc_histos_%d_rad.root", run);                      
  TString data_filename = Form("heep_data_histos_%d_combined.root",run); 

  TString simc_f = root_dir + simc_filename;
  TString data_f = root_dir + data_filename;
  
  //Where to store plots
  string plots_dir = "./";
  string plots_path;
  
  //Open SIMC/data ROOT files;
  TFile *simc_file = new TFile(simc_f);
  TFile *data_file = new TFile(data_f);

  //---------------Target ----------------
  //Define SIMC histos ('h'-->hadron arm,  'e'-->electron arm)
  
  TH1F *simc_xtar =  0;
  TH1F *simc_ytarH =  0;
  TH1F *simc_ztarH =  0;

  TH1F *simc_ytarP =  0;
  TH1F *simc_ztarP =  0;  

  //Define data histos
  TH1F *data_xtarH = 0;
  TH1F *data_ytarH = 0;
  TH1F *data_ztarH = 0;

  TH1F *data_xtarP = 0;                                                                                                     
  TH1F *data_ytarP = 0;                                                                                                                                   
  TH1F *data_ztarP = 0; 

  //---------------Target Reconstruction Variables----------------
  //Define SIMC histos ('h'-->hadron arm,  'e'-->electron arm)
  TH1F *simc_eytar =  0;
  TH1F *simc_exptar =  0;
  TH1F *simc_eyptar =  0;
  TH1F *simc_edelta =  0;

  TH1F *simc_hytar =  0;
  TH1F *simc_hxptar =  0;
  TH1F *simc_hyptar =  0;
  TH1F *simc_hdelta =  0;

  //Define data histos
  TH1F *data_eytar = 0;
  TH1F *data_exptar =  0;
  TH1F *data_eyptar =  0;
  TH1F *data_edelta =  0;

  TH1F *data_hytar = 0;
  TH1F *data_hxptar =  0;
  TH1F *data_hyptar =  0;
  TH1F *data_hdelta =  0;

  //-----------------------------------------------------------
 
  //--------------FOCAL PLANE VARIABLES------------------------

 //Define SIMC histos ('h'-->hadron arm,  'e'-->electron arm)
  TH1F *simc_exfp =  0;
  TH1F *simc_eyfp =  0;
  TH1F *simc_expfp =  0;
  TH1F *simc_eypfp =  0;

  TH1F *simc_hxfp =  0;
  TH1F *simc_hyfp =  0;
  TH1F *simc_hxpfp =  0;
  TH1F *simc_hypfp =  0;
  
  //Define data histos
  TH1F *data_exfp =  0;
  TH1F *data_eyfp =  0;
  TH1F *data_expfp =  0;
  TH1F *data_eypfp =  0;

  TH1F *data_hxfp =  0;
  TH1F *data_hyfp =  0;
  TH1F *data_hxpfp =  0;
  TH1F *data_hypfp =  0;

  //--------------------------------------------------------------

  //-------------------------KINEMATICS---------------------------
  TH1F *simc_Q2 =  0;
  TH1F *simc_omega =  0;
  TH1F *simc_W =  0;
  TH1F *simc_thq = 0;

  TH1F *simc_xbj = 0;
  TH1F *simc_th_elec = 0;                                  
  TH1F *simc_kf = 0;  
  TH1F *simc_emiss = 0;

  //Kinematics 2
  TH1F *simc_Pm = 0;
  TH1F *simc_Pf = 0;
  TH1F *simc_th_prot = 0;
  TH1F *simc_q = 0;    //q-vector magnitude
  TH1F *simc_thpq = 0;
  TH1F *simc_Pmx = 0;
  TH1F *simc_Pmy = 0;
  TH1F *simc_Pmz = 0;

  //Define data histos
  TH1F *data_Q2 =  0;
  TH1F *data_omega =  0;
  TH1F *data_W =  0;
  TH1F *data_thq = 0;

  TH1F *data_xbj = 0;
  TH1F *data_th_elec = 0;
  TH1F *data_kf = 0;
  TH1F *data_emiss = 0;

   //Kinematics 2
  TH1F *data_Pm = 0;
  TH1F *data_Pf = 0;
  TH1F *data_th_prot = 0;
  TH1F *data_q = 0;    //q-vector magnitude
  TH1F *data_thpq = 0;
  TH1F *data_Pmx = 0;
  TH1F *data_Pmy = 0;
  TH1F *data_Pmz = 0;

  //-----------------------
  TH1F *simc_ztar_diff = 0; 
  TH2F *simc_HMS_Coll = 0;
  
  TH1F *data_ztar_diff = 0; 
  TH2F *data_HMS_Coll = 0;

  TH1F *data_CoinTime = 0;
  TH1F *data_pid_eCal = 0;

  //---------------------------------------------------------------

 //change to simc_file
  simc_file->cd();

  simc_file->GetObject("H_ztar_diff_sys", simc_ztar_diff);
  simc_file->GetObject("H_hXColl_vs_hYColl_sys", simc_HMS_Coll);

  
  
  //----------Get Target Histograms------------------
  //Get Histogram objects from SIMC rootfile
  simc_file->GetObject("H_hx_tar", simc_xtar);

  simc_file->GetObject("H_hy_tar", simc_ytarH);
  simc_file->GetObject("H_hz_tar", simc_ztarH);

  simc_file->GetObject("H_ey_tar", simc_ytarP);  
  simc_file->GetObject("H_ez_tar", simc_ztarP); 

  
  //change to data_file
  data_file->cd();

  data_file->GetObject("H_ztar_diff_sys", data_ztar_diff);
  data_file->GetObject("H_hXColl_vs_hYColl_sys", data_HMS_Coll); 
  data_file->GetObject("H_pcal_etotTrkNorm_sys", data_pid_eCal); 
  data_file->GetObject("H_ctime_sys", data_CoinTime);
  
  //Get Histogram objects from data rootfile
  data_file->GetObject("H_hx_tar", data_xtarH);
  data_file->GetObject("H_hy_tar", data_ytarH);
  data_file->GetObject("H_hz_tar", data_ztarH);

  data_file->GetObject("H_ex_tar", data_xtarP); 
  data_file->GetObject("H_ey_tar", data_ytarP);                        
  data_file->GetObject("H_ez_tar", data_ztarP); 
  


  //-----------------------------------------------------------------


  //---------------------------------------------------------------

 //change to simc_file
  simc_file->cd();

  //----------Get Target Reconstructed Histograms------------------
  //Get Histogram objects from SIMC rootfile
  simc_file->GetObject("H_eytar", simc_eytar);
  simc_file->GetObject("H_exptar", simc_exptar);
  simc_file->GetObject("H_eyptar", simc_eyptar);
  simc_file->GetObject("H_edelta_sys", simc_edelta);

  simc_file->GetObject("H_hytar", simc_hytar);
  simc_file->GetObject("H_hxptar", simc_hxptar);
  simc_file->GetObject("H_hyptar", simc_hyptar);
  simc_file->GetObject("H_hdelta_sys", simc_hdelta);

  simc_file->GetObject("H_Pmx_Lab", simc_Pmx);
  simc_file->GetObject("H_Pmy_Lab", simc_Pmy);
  simc_file->GetObject("H_Pmz_Lab", simc_Pmz);



  //change to data_file
  data_file->cd();

  //Get Histogram objects from data rootfile
  data_file->GetObject("H_eytar", data_eytar);
  data_file->GetObject("H_exptar", data_exptar);
  data_file->GetObject("H_eyptar", data_eyptar);
  data_file->GetObject("H_edelta_sys", data_edelta);
  
  data_file->GetObject("H_hytar", data_hytar);
  data_file->GetObject("H_hxptar", data_hxptar);
  data_file->GetObject("H_hyptar", data_hyptar);
  data_file->GetObject("H_hdelta_sys", data_hdelta);


  //-----------------------------------------------------------------

  

  //---------------Get FOCAL PLANE Histograms------------------------

   //change to simc_file
  simc_file->cd();

  //Get Histogram objects from SIMC rootfile
  simc_file->GetObject("H_exfp", simc_exfp);
  simc_file->GetObject("H_eyfp", simc_eyfp);
  simc_file->GetObject("H_expfp", simc_expfp);
  simc_file->GetObject("H_eypfp", simc_eypfp);

  simc_file->GetObject("H_hxfp", simc_hxfp);
  simc_file->GetObject("H_hyfp", simc_hyfp);
  simc_file->GetObject("H_hxpfp", simc_hxpfp);
  simc_file->GetObject("H_hypfp", simc_hypfp);

  //change to data_file
  data_file->cd();

  //Get Histogram objects from data rootfile
  data_file->GetObject("H_exfp", data_exfp);
  data_file->GetObject("H_eyfp", data_eyfp);
  data_file->GetObject("H_expfp", data_expfp);
  data_file->GetObject("H_eypfp", data_eypfp);

  data_file->GetObject("H_hxfp", data_hxfp);
  data_file->GetObject("H_hyfp", data_hyfp);
  data_file->GetObject("H_hxpfp", data_hxpfp);
  data_file->GetObject("H_hypfp", data_hypfp);

  //--------------------------------------------------------------
  
  //------------------Get KINEMATICS VARIABLES--------------------

   //change to simc_file
  simc_file->cd();

  //Get Histogram objects from SIMC rootfile
  simc_file->GetObject("H_Q2_sys", simc_Q2);
  simc_file->GetObject("H_omega", simc_omega);
  simc_file->GetObject("H_W", simc_W);
  simc_file->GetObject("H_theta_q", simc_thq);

  simc_file->GetObject("H_xbj", simc_xbj);
  simc_file->GetObject("H_theta_elec", simc_th_elec);
  simc_file->GetObject("H_kf", simc_kf);
  simc_file->GetObject("H_Em_sys", simc_emiss);

  simc_file->GetObject("H_Pm", simc_Pm);
  simc_file->GetObject("H_Pf", simc_Pf);
  simc_file->GetObject("H_theta_prot", simc_th_prot);
  simc_file->GetObject("H_q", simc_q);
  simc_file->GetObject("H_theta_pq", simc_thpq);


  //change to data_file
  data_file->cd();
  
  //Get Histogram objects from data rootfile
  data_file->GetObject("H_Q2_sys", data_Q2);
  data_file->GetObject("H_omega", data_omega);
  data_file->GetObject("H_W", data_W);
  data_file->GetObject("H_theta_q", data_thq);

  
  data_file->GetObject("H_xbj", data_xbj);
  data_file->GetObject("H_theta_elec", data_th_elec);
  data_file->GetObject("H_kf", data_kf);
  data_file->GetObject("H_Em_sys", data_emiss);

  data_file->GetObject("H_Pm", data_Pm);
  data_file->GetObject("H_Pf", data_Pf);
  data_file->GetObject("H_theta_prot", data_th_prot);
  data_file->GetObject("H_q", data_q);
  data_file->GetObject("H_theta_pq", data_thpq);

  data_file->GetObject("H_Pmx_Lab", data_Pmx);
  data_file->GetObject("H_Pmy_Lab", data_Pmy);
  data_file->GetObject("H_Pmz_Lab", data_Pmz);


  //Plot the SHMS Recon Ratio
  //hist_ratio(data_eytar, simc_eytar, "SHMS Y_{tar} [cm]", "Counts / mC", "SHMS Y_{tar}");
  //hist_ratio(data_exptar, simc_exptar, "SHMS X'_{tar} [rad]", "Counts / mC", "SHMS X'_{tar}");
  //hist_ratio(data_eyptar, simc_eyptar, "SHMS Y'_{tar} [rad]", "Counts / mC", "SHMS Y'_{tar}");
  //hist_ratio(data_edelta, simc_edelta, "SHMS #delta [%]", "Counts / mC", "SHMS #delta");

  //Plot the HMS Recon Ratio
  //hist_ratio(data_hytar, simc_hytar, "HMS Y_{tar} [cm]", "Counts / mC", "HMS Y_{tar}");
  //hist_ratio(data_hxptar, simc_hxptar, "HMS X'_{tar} [rad]", "Counts / mC", "HMS X'_{tar}");
  //hist_ratio(data_hyptar, simc_hyptar, "HMS Y'_{tar} [rad]", "Counts / mC", "HMS Y'_{tar}");
  //hist_ratio(data_hdelta, simc_hdelta, "HMS #delta [%]", "Counts / mC", "HMS #delta");
  
  
  //Plot Kinematic Variables
  //hist_ratio(data_W, simc_W, "Invariant Mass, W [GeV]", "Counts", "Invariant Mass");

  
  //ONLY COMPARE HISTOS  
   //Plot the SHMS Recon
  
  //compare_hist(data_eytar, simc_eytar, "SHMS Y_{tar} [cm]", "Counts", "SHMS Y_{tar}");
  //compare_hist(data_exptar, simc_exptar, "SHMS X'_{tar} [rad]", "Counts", "SHMS X'_{tar}");
  //compare_hist(data_eyptar, simc_eyptar, "SHMS Y'_{tar} [rad]", "Counts", "SHMS Y'_{tar}");
  //compare_hist(data_edelta, simc_edelta, "SHMS #delta [%]", "Counts", "SHMS #delta");

  //Plot the HMS Recon
  //compare_hist(data_hytar, simc_hytar, "HMS Y_{tar} [cm]", "Counts", "HMS Y_{tar}");
  //compare_hist(data_hxptar, simc_hxptar, "HMS X'_{tar} [rad]", "Counts", "HMS X'_{tar}");
  //compare_hist(data_hyptar, simc_hyptar, "HMS Y'_{tar} [rad]", "Counts", "HMS Y'_{tar}");
  //compare_hist(data_hdelta, simc_hdelta, "HMS #delta [%]", "Counts", "HMS #delta");
  

  //------Plot Event Selection Cuts----
  compare_hist(data_W, simc_W, "Invariant Mass W [GeV]", "Counts / mC", Form("H(e,e'p) Elastic Run %d: Invariant Mass", run));
  
  //Also plot collimator
  //plot_hist(data_pid_eCal, "SHMS Calorimeter E_{dep}/P_{trk}", "Counts / mC", "SHMS Calorimeter Total Normalized Energy", "logy");
  //plot_hist(data_CoinTime, "Coincidence Time [ns]", "Counts / mC", "Coincidence Time", "logy");

  //ONLY PLOT DATA-2-SIMC COMPARISONS (NOT RATIOS)
  //compare_hist(data_hdelta, simc_hdelta, "HMS #delta [%]", "Counts / mC", "HMS #delta");
  //compare_hist(data_edelta, simc_edelta, "SHMS #delta [%]", "Counts / mC", "SHMS #delta");

  compare_hist(data_Q2, simc_Q2, "Q^{2} [GeV^{2}]", "Counts / mC", Form("H(e,e'p) Elastic Run %d: 4-Momentum Transfer, Q^{2}", run));
  //compare_hist(data_thnq, simc_thnq_fsi, "#theta_{nq} [deg]", "Charge Normalized Counts", "Neutron Recoil Angles, #theta_{nq}");
  //compare_hist(data_ztar_diff, simc_ztar_diff, "Z_{tar} Difference [cm]", "Counts / mC", "Z_{tar} Difference", "logy");
  //compare_hist(data_emiss, simc_emiss, "Missing Energy, E_{m} [GeV]", "Counts / mC", "Missing Energy");


  //------Plott elecntron angle/momentum
  compare_hist(data_th_elec, simc_th_elec, "#theta_{e} [deg]", "Charge Normalized Counts", "Electron Scattered Angles, #theta_{e}");
  compare_hist(data_kf, simc_kf, "k_{f} [GeV/c]", "Charge Normalized Counts", "Electron Scattered Momentum, k_{f}");
  
}
