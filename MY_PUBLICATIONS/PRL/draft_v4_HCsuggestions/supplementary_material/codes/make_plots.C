void make_plot(int set, int pm, string model, string rad)
{

  //gROOT->SetBatch(kTRUE);  
  gStyle->SetOptStat(1001111);

  string thnq_set = "all_thnq_Q2_4to5GeV";
  //string thnq_set = "thnq35_Q2_4to5GeV";
  //string thnq_set = "thnq45_Q2_4to5GeV";
  //string thnq_set = "thnq75_Q2_4to5GeV";
  
  
  TString root_dir = Form("../root_files/pm%d_fsiXsec_set%d_%s/", pm, set, thnq_set.c_str());
  TString simc_fname = Form("deep_simc_histos_pm%d_lagetfsi_rad_set%d.root", pm, set);
  TString data_fname = Form("deep_data_histos_pm%d_set%d_combined.root", pm, set);
    
  TString simc_filename_fsi = root_dir+simc_fname;
  TString data_filename = root_dir+data_fname;
  
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

  //Set SIMC Histo Aesthetics          
  simc_ztar_diff_fsi->SetFillColorAlpha(kRed, 0.35);                                                                                                                                   
  simc_ztar_diff_fsi->SetFillStyle(3004);                                                                                                                                            
  simc_ztar_diff_fsi->SetLineColor(kRed); 


  //------SIMC Spectrometer Acceptance------
  simc_file_fsi->GetObject("H_eytar", simc_eytar_fsi);
  simc_file_fsi->GetObject("H_exptar", simc_exptar_fsi);
  simc_file_fsi->GetObject("H_eyptar", simc_eyptar_fsi);
  simc_file_fsi->GetObject("H_edelta_sys", simc_edelta_fsi);

  simc_file_fsi->GetObject("H_hytar", simc_hytar_fsi);
  simc_file_fsi->GetObject("H_hxptar", simc_hxptar_fsi);
  simc_file_fsi->GetObject("H_hyptar", simc_hyptar_fsi);
  simc_file_fsi->GetObject("H_hdelta_sys", simc_hdelta_fsi);

  //Set SIMC Histo Aesthetics
  simc_eytar_fsi->SetFillColorAlpha(kRed, 0.35);
  simc_eytar_fsi->SetFillStyle(3004);
  simc_eytar_fsi->SetLineColor(kRed);
  
  simc_exptar_fsi->SetFillColorAlpha(kRed, 0.35);
  simc_exptar_fsi->SetFillStyle(3004);
  simc_exptar_fsi->SetLineColor(kRed);
    
  simc_eyptar_fsi->SetFillColorAlpha(kRed, 0.35);
  simc_eyptar_fsi->SetFillStyle(3004);
  simc_eyptar_fsi->SetLineColor(kRed);

  simc_edelta_fsi->SetFillColorAlpha(kRed, 0.35);
  simc_edelta_fsi->SetFillStyle(3004);
  simc_edelta_fsi->SetLineColor(kRed);

  simc_hytar_fsi->SetFillColorAlpha(kRed, 0.35);
  simc_hytar_fsi->SetFillStyle(3004);
  simc_hytar_fsi->SetLineColor(kRed);

  simc_hxptar_fsi->SetFillColorAlpha(kRed, 0.35);
  simc_hxptar_fsi->SetFillStyle(3004);
  simc_hxptar_fsi->SetLineColor(kRed);

  simc_hyptar_fsi->SetFillColorAlpha(kRed, 0.35);
  simc_hyptar_fsi->SetFillStyle(3004);
  simc_hyptar_fsi->SetLineColor(kRed);

  simc_hdelta_fsi->SetFillColorAlpha(kRed, 0.35);
  simc_hdelta_fsi->SetFillStyle(3004);
  simc_hdelta_fsi->SetLineColor(kRed);


  //------SIMC ANALYSIS PLOTS-----
  simc_file_fsi->GetObject("H_Q2_sys", simc_Q2_fsi);
  simc_file_fsi->GetObject("H_theta_nq_sys", simc_thnq_fsi);
  simc_file_fsi->GetObject("H_Em_nuc_sys", simc_emiss_fsi);

  simc_file_fsi->GetObject("H_Pm", simc_Pm_fsi);

  //Set SIMC Histo Aesthetics
  simc_Q2_fsi->SetFillColorAlpha(kRed, 0.35);
  simc_Q2_fsi->SetFillStyle(3004);
  simc_Q2_fsi->SetLineColor(kRed);

  simc_thnq_fsi->SetFillColorAlpha(kRed, 0.35);
  simc_thnq_fsi->SetFillStyle(3004);
  simc_thnq_fsi->SetLineColor(kRed);
  
  simc_emiss_fsi->SetFillColorAlpha(kRed, 0.35);
  simc_emiss_fsi->SetFillStyle(3004);
  simc_emiss_fsi->SetLineColor(kRed);

  simc_Pm_fsi->SetFillColorAlpha(kRed, 0.35);
  simc_Pm_fsi->SetFillStyle(3004);
  simc_Pm_fsi->SetLineColor(kRed);

  //Get Histogram objects from DATA rootfile
  //change to data_file
  data_file->cd();
  
  data_file->GetObject("H_ztar_diff_sys", data_ztar_diff);
  data_file->GetObject("H_hXColl_vs_hYColl_sys", data_HMS_Coll); 
  data_file->GetObject("H_pcal_etotTrkNorm_sys", data_pid_eCal); 
  data_file->GetObject("H_ctime_sys", data_CoinTime); 

  //Set DATA Histo Aesthetics                                                                                                                         
  data_ztar_diff->SetLineColor(kBlue);  
  data_ztar_diff->SetLineWidth(2); 
    
  data_pid_eCal->SetLineColor(kBlue); 
  data_pid_eCal->SetLineWidth(2); 
  data_pid_eCal->SetFillColorAlpha(kBlue, 0.35);
  data_pid_eCal->SetFillStyle(3004);
  
  data_CoinTime->SetLineColor(kBlue); 
  data_CoinTime->SetLineWidth(2); 
  data_CoinTime->SetFillColorAlpha(kBlue, 0.35);
  data_CoinTime->SetFillStyle(3004);
  
  //-------DATA Spectrometer Acceptance-----
  data_file->GetObject("H_eytar", data_eytar);
  data_file->GetObject("H_exptar", data_exptar);
  data_file->GetObject("H_eyptar", data_eyptar);
  data_file->GetObject("H_edelta_sys", data_edelta);
  
  data_file->GetObject("H_hytar", data_hytar);
  data_file->GetObject("H_hxptar", data_hxptar);
  data_file->GetObject("H_hyptar", data_hyptar);
  data_file->GetObject("H_hdelta_sys", data_hdelta);

  //Set data Histo Aesthetics
  data_eytar->SetLineColor(kBlue);
  data_eytar->SetLineWidth(2);
  data_exptar->SetLineColor(kBlue);
  data_exptar->SetLineWidth(2);
  data_eyptar->SetLineColor(kBlue);
  data_eyptar->SetLineWidth(2);
  data_edelta->SetLineColor(kBlue);
  data_edelta->SetLineWidth(2);

  data_hytar->SetLineColor(kBlue);
  data_hytar->SetLineWidth(2);
  data_hxptar->SetLineColor(kBlue);
  data_hxptar->SetLineWidth(2);
  data_hyptar->SetLineColor(kBlue);
  data_hyptar->SetLineWidth(2);
  data_hdelta->SetLineColor(kBlue);
  data_hdelta->SetLineWidth(2);

  //-----Get Data Analysis Histograms
  data_file->GetObject("H_Q2_sys", data_Q2);
  data_file->GetObject("H_theta_nq_sys", data_thnq);
  data_file->GetObject("H_Em_nuc_sys", data_emiss);
  
  data_file->GetObject("H_Pm", data_Pm);

  
  //Set data Histo Aesthetics
  data_Q2->SetLineColor(kBlue);
  data_Q2->SetLineWidth(2);

  data_thnq->SetLineColor(kBlue);
  data_thnq->SetLineWidth(2);
  
  data_emiss->SetLineColor(kBlue);
  data_emiss->SetLineWidth(2);

  data_Pm->SetLineColor(kBlue);
  data_Pm->SetLineWidth(2);


  //====Plot Spectrometer Acceptance
  TCanvas *c1 = new TCanvas("c1", "Electron Arm: Target Reconstruction", 5000, 3000);
  c1->Divide(2,2);
  
  c1->cd(1);
  data_eytar->Draw("samesE0");
  simc_eytar_fsi->Draw("sameshistE0");
  
  c1->cd(2);
  data_exptar->Draw("samesE0");
  simc_exptar_fsi->Draw("sameshistE0");
  
  c1->cd(3);
  data_eyptar->Draw("samesE0");
  simc_eyptar_fsi->Draw("sameshistE0");
  
  c1->cd(4);
  data_edelta->Draw("samesE0");
  simc_edelta_fsi->Draw("sameshistE0");
  

  TCanvas *c2 = new TCanvas("c2", "Hadron Arm: Target Reconstruction", 5000, 3000);
  c2->Divide(2,2);
  
  c2->cd(1);
  data_hytar->Draw("samesE0");
  simc_hytar_fsi->Draw("sameshistE0");
  
  c2->cd(2);
  data_hxptar->Draw("samesE0");
  simc_hxptar_fsi->Draw("sameshistE0");
  
  c2->cd(3);
  data_hyptar->Draw("samesE0");
  simc_hyptar_fsi->Draw("sameshistE0");
  
  c2->cd(4);
  data_hdelta->Draw("samesE0");
  simc_hdelta_fsi->Draw("sameshistE0");  

  //-------Plot Analysis Cuts----
  
  TCanvas *c3 = new TCanvas("c3", "Analysis Cuts", 5000, 4000);
  c3->Divide(3,2);  //Em, ztar_diff, Q2, thnq, coin_tim, shms_cal

  c3->cd(1);
  data_emiss->Draw("samesE0");
  simc_emiss_fsi->Draw("sameshistE0");  

  c3->cd(2);
  data_ztar_diff->Draw("samesE0");
  simc_ztar_diff_fsi->Draw("sameshistE0");  

  c3->cd(3);
  data_Q2->Draw("samesE0");
  simc_Q2_fsi->Draw("sameshistE0");  

  c3->cd(4);
  data_thnq->Draw("samesE0");
  simc_thnq_fsi->Draw("sameshistE0");  

  c3->cd(5);
  data_CoinTime->Draw("sameshistE0");

  c3->cd(6);
  data_pid_eCal->Draw("sameshistE0");

  //-----Plot HMS Collimator
  TCanvas *c4 = new TCanvas("c4", "Analysis Cuts Coll.", 2000, 4000);
  c4->Divide(2,1);  //Em, ztar_diff, Q2, thnq, coin_tim, shms_cal

  c4->cd(1);
  data_HMS_Coll->Draw("colz");

  c4->cd(2);
  simc_HMS_Coll_fsi->Draw("colz");

}
