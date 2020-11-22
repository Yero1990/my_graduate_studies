#include "../plotting_utilities/histUtils.C"   //load utilities code for histo plotting

void plot_combined_sets(int thnq=0)
{
 gStyle->SetOptStat(0000000);

 //set filenames for plotting kinematics
 /*
 string fname_80 = Form("root_files_combined/%d_deg/deep_data_histos_pm80_set1_combined.root", thnq);
 string fname_580_set1 = Form("root_files_combined/%d_deg/deep_data_histos_pm580_set1_combined.root",thnq);
 string fname_580_set2 = Form("root_files_combined/%d_deg/deep_data_histos_pm580_set2_combined.root", thnq);

 string fname_750_set1 = Form("root_files_combined/%d_deg/deep_data_histos_pm750_set1_combined.root", thnq);
 string fname_750_set2 = Form("root_files_combined/%d_deg/deep_data_histos_pm750_set2_combined.root", thnq);
 string fname_750_set3 = Form("root_files_combined/%d_deg/deep_data_histos_pm750_set3_combined.root", thnq);
 */

 /*
 //set filenames for plotting Missing Momentum yield (rad. corrected)
 string fname_80 = Form("../../root_files/pm80_fsiXsec_set1_thnq%d_Q2_4to5/deep_data_histos_pm80_set1_combined_radcorr.root", thnq);

 string fname_580_set1 = Form("../../root_files/pm580_fsiXsec_set1_thnq%d_Q2_4to5/deep_data_histos_pm580_set1_combined_radcorr.root",thnq);
 string fname_580_set2 = Form("../../root_files/pm580_fsiXsec_set2_thnq%d_Q2_4to5/deep_data_histos_pm580_set2_combined_radcorr.root",thnq);

 string fname_750_set1 = Form("../../root_files/pm750_fsiXsec_set1_thnq%d_Q2_4to5/deep_data_histos_pm750_set1_combined_radcorr.root",thnq);
 string fname_750_set2 = Form("../../root_files/pm750_fsiXsec_set2_thnq%d_Q2_4to5/deep_data_histos_pm750_set2_combined_radcorr.root",thnq);
 string fname_750_set3 = Form("../../root_files/pm750_fsiXsec_set3_thnq%d_Q2_4to5/deep_data_histos_pm750_set3_combined_radcorr.root",thnq);
 */

 
 //set filenames for plotting Missing Momentum yield (rad. corrected)
 string fname_80 = Form("../../root_files/pm80_fsiXsec_set1_thnq%d_Q2_4to5/deep_simc_histos_pm80_lagetfsi_norad_set1.root", thnq);

 string fname_580_set1 = Form("../../root_files/pm580_fsiXsec_set1_thnq%d_Q2_4to5/deep_simc_histos_pm580_lagetfsi_norad_set1.root",thnq);
 string fname_580_set2 = Form("../../root_files/pm580_fsiXsec_set2_thnq%d_Q2_4to5/deep_simc_histos_pm580_lagetfsi_norad_set2.root",thnq);

 string fname_750_set1 = Form("../../root_files/pm750_fsiXsec_set1_thnq%d_Q2_4to5/deep_simc_histos_pm750_lagetfsi_norad_set1.root",thnq);
 string fname_750_set2 = Form("../../root_files/pm750_fsiXsec_set2_thnq%d_Q2_4to5/deep_simc_histos_pm750_lagetfsi_norad_set2.root",thnq);
 string fname_750_set3 = Form("../../root_files/pm750_fsiXsec_set3_thnq%d_Q2_4to5/deep_simc_histos_pm750_lagetfsi_norad_set3.root",thnq);
 
 
 TFile *f_80 = new TFile(fname_80.c_str());
 TFile *f_580_1 = new TFile(fname_580_set1.c_str());
 TFile *f_580_2 = new TFile(fname_580_set2.c_str());

 TFile *f_750_1 = new TFile(fname_750_set1.c_str());
 TFile *f_750_2 = new TFile(fname_750_set2.c_str());
 TFile *f_750_3 = new TFile(fname_750_set3.c_str());

 //Create Hist object
 TH1F *hist_80 = 0;
 TH1F *hist_580_set1 = 0;
 TH1F *hist_580_set2 = 0;
 TH1F *hist_750_set1 = 0;
 TH1F *hist_750_set2 = 0;
 TH1F *hist_750_set3 = 0;
 
 TString hist_name, hist_title, hist_xlabel;
 //hist_name = "H_Q2_sys", hist_title = "4-Momentum Transfer, Q^{2}", hist_xlabel = "Q^{2} [(GeV/c)^{2}]";
 //hist_name = "H_xbj", hist_title = "x-Bjorken", hist_xlabel = "#it{x}_{Bj}";
 //hist_name = "H_omega", hist_title = "Energy Transfer, #omega", hist_xlabel = "#omega [GeV]";
 //hist_name = "H_Pf", hist_title = "Final Proton Momentum, #it{p}_{f}", hist_xlabel = "#it{p}_{f} [GeV/c]";
 hist_name = "H_Pm_ps", hist_title = "Missing Momentum Phase Space", hist_xlabel = "#it{p}_{r} [GeV/c]";
 
 //GET HISTOS FROM EACH SET
 f_80->cd();
 f_80->GetObject(hist_name.Data(), hist_80);
 
 f_580_1->cd();
 f_580_1->GetObject(hist_name.Data(), hist_580_set1);

 f_580_2->cd();
 f_580_2->GetObject(hist_name.Data(), hist_580_set2);
 
 f_750_1->cd();
 f_750_1->GetObject(hist_name.Data(), hist_750_set1);
  
 f_750_2->cd();
 f_750_2->GetObject(hist_name.Data(), hist_750_set2);

 f_750_3->cd();
 f_750_3->GetObject(hist_name.Data(), hist_750_set3);

 combine_sets(hist_80, hist_580_set1, hist_580_set2, hist_750_set1, hist_750_set2, hist_750_set3, hist_xlabel, "arb. units", hist_title, 1., 0);
 //combine_sets(hist_80, 0, 0, 0, 0, 0, "arb. units", hist_title, 500.);
 
}
