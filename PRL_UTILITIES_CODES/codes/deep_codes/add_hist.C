#include "../plotting_utilities/histUtils.C"   //load utilities code for histo plotting

//Script to make comparison between SIMC and Commissioning Data from HallC Spring 2018
//The comparisons are all charge normalized and corrected for all inefficiencies

void add_hist(int thnq)
{

  //gROOT->SetBatch(kTRUE);  
  //gStyle->SetOptStat(1001111);
  gStyle->SetOptStat(0000000);
  string thnq_set;
  string root_dir_1;
  string root_dir_2;
  string root_dir_3;
  string root_dir_4;
  string root_dir_5;


  if(thnq==-1){
    //thnq_set = "all_thnq_Q2_4to5";
    thnq_set = "all_thnq_Q2_4to5_WideEmBINS";
  }
  else{thnq_set = Form("thnq%d_Q2_4to5", thnq);}

  
  
  if(thnq==-1){
    root_dir_1 = Form("../../root_files/pm580_fsiXsec_set1_%s/deep_data_histos_pm580_set1_combined.root", thnq_set.c_str());
    root_dir_2 = Form("../../root_files/pm580_fsiXsec_set2_%s/deep_data_histos_pm580_set2_combined.root", thnq_set.c_str());
    root_dir_3 = Form("../../root_files/pm750_fsiXsec_set1_%s/deep_data_histos_pm750_set1_combined.root", thnq_set.c_str());
    root_dir_4 = Form("../../root_files/pm750_fsiXsec_set2_%s/deep_data_histos_pm750_set2_combined.root", thnq_set.c_str());
    root_dir_5 = Form("../../root_files/pm750_fsiXsec_set3_%s/deep_data_histos_pm750_set3_combined.root", thnq_set.c_str());
  }

  
  TFile *data_file_580_set1 = new TFile(root_dir_1.c_str());
  TFile *data_file_580_set2 = new TFile(root_dir_2.c_str());

  TFile *data_file_750_set1 = new TFile(root_dir_3.c_str());
  TFile *data_file_750_set2 = new TFile(root_dir_4.c_str());
  TFile *data_file_750_set3 = new TFile(root_dir_5.c_str());
  

  TH1F *data_Em_580_set1 = 0;
  TH1F *data_Em_580_set2 = 0;
  TH1F *data_Em_750_set1 = 0;
  TH1F *data_Em_750_set2 = 0;
  TH1F *data_Em_750_set3 = 0;
  


  data_file_580_set1->cd();
  data_file_580_set1->GetObject("H_Em_nuc_sys", data_Em_580_set1);  

  data_file_580_set2->cd();
  data_file_580_set2->GetObject("H_Em_nuc_sys", data_Em_580_set2);  

  data_file_750_set1->cd();
  data_file_750_set1->GetObject("H_Em_nuc_sys", data_Em_750_set1);  

  data_file_750_set2->cd();
  data_file_750_set2->GetObject("H_Em_nuc_sys", data_Em_750_set2);  
 
  data_file_750_set3->cd();
  data_file_750_set3->GetObject("H_Em_nuc_sys", data_Em_750_set3);  
   

  //============================================================================================

  /*
  TCanvas *c1 = new TCanvas("c1", "", 500, 500);
  data_Em_580_set1->Add(data_Em_580_set2);
  data_Em_580_set1->Draw();

  Double_t err_580_dummy = 0.0;
  Double_t err_580_data = 0.0;

  Double_t I_580_dummy = data_Em_580_set1->IntegralAndError(data_Em_580_set1->FindBin(-0.06), data_Em_580_set1->FindBin(-0.03), err_580_dummy);
  Double_t I_580_data = data_Em_580_set1->IntegralAndError(data_Em_580_set1->FindBin(-0.01), data_Em_580_set1->FindBin(0.02), err_580_data);

  cout << "Al_dummy yield = " << I_580_dummy << " +/- " << err_580_dummy << endl;
  cout << "deut+dummy yield = " << I_580_data << " +/- " << err_580_data << endl;
  cout << "fractional dummy yield / (data+dummy) = " <<  (I_580_dummy /  I_580_data) * 100 << endl; 
  */
  //TCanvas *c2 = new TCanvas("c2", "", 500, 500);
  data_Em_750_set1->Add(data_Em_750_set2);
  data_Em_750_set1->Add(data_Em_750_set3);
  data_Em_750_set1->Draw();
  
  Double_t err_750_dummy = 0.0;
  Double_t err_750_data = 0.0;

  Double_t I_750_dummy = data_Em_750_set1->IntegralAndError(data_Em_750_set1->FindBin(-0.06), data_Em_750_set1->FindBin(-0.03), err_750_dummy);
  Double_t I_750_data = data_Em_750_set1->IntegralAndError(data_Em_750_set1->FindBin(-0.01), data_Em_750_set1->FindBin(0.02), err_750_data);

  cout << "Al_dummy yield = " << I_750_dummy << " +/- " << err_750_dummy << endl;
  cout << "deut+dummy yield = " << I_750_data << " +/- " << err_750_data << endl;
  cout << "fractional dummy yield / (data+dummy) = " <<  (I_750_dummy /  I_750_data) * 100 << endl; 
  
}
