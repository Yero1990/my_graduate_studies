#include "../plotting_utilities/histUtils.C"   //load utilities code for histo plotting

//Script to make comparison between SIMC and Commissioning Data from HallC Spring 2018
//The comparisons are all charge normalized and corrected for all inefficiencies

void add_hist(int thnq)
{

  //gROOT->SetBatch(kTRUE);  
  //gStyle->SetOptStat(1001111);
  gStyle->SetOptStat(0000000);
  string thnq_set;
  string root_dir_0;
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
    //root_dir_0 = Form("../../root_files/pm80_fsiXsec_set1_%s/deep_data_histos_pm80_set1_combined.root", thnq_set.c_str());
    //root_dir_1 = Form("../../root_files/pm580_fsiXsec_set1_%s/deep_data_histos_pm580_set1_combined.root", thnq_set.c_str());
    //root_dir_2 = Form("../../root_files/pm580_fsiXsec_set2_%s/deep_data_histos_pm580_set2_combined.root", thnq_set.c_str());
    //root_dir_3 = Form("../../root_files/pm750_fsiXsec_set1_%s/deep_data_histos_pm750_set1_combined.root", thnq_set.c_str());
    //root_dir_4 = Form("../../root_files/pm750_fsiXsec_set2_%s/deep_data_histos_pm750_set2_combined.root", thnq_set.c_str());
    //root_dir_5 = Form("../../root_files/pm750_fsiXsec_set3_%s/deep_data_histos_pm750_set3_combined.root", thnq_set.c_str());

    root_dir_1 = "../../root_files/deep_data_histos_pm580_set1_combined.root";
    root_dir_2 = "../../root_files/deep_data_histos_pm580_set2_combined.root";
    
  }

  //TFile *data_file_80_set1 = new TFile(root_dir_0.c_str());

  TFile *data_file_580_set1 = new TFile(root_dir_1.c_str());
  TFile *data_file_580_set2 = new TFile(root_dir_2.c_str());

  //TFile *data_file_750_set1 = new TFile(root_dir_3.c_str());
  //TFile *data_file_750_set2 = new TFile(root_dir_4.c_str());
  //TFile *data_file_750_set3 = new TFile(root_dir_5.c_str());
  
  TH1F *data_Em_80_set1 = 0;
  TH1F *data_Em_580_set1 = 0;
  TH1F *data_Em_580_set2 = 0;
  TH1F *data_Em_750_set1 = 0;
  TH1F *data_Em_750_set2 = 0;
  TH1F *data_Em_750_set3 = 0;
  TH1F *data_eztar_580_set1 = 0;
  TH1F *data_eztar_580_set2 = 0;
 
  //data_file_80_set1->cd();
  //data_file_80_set1->GetObject("H_Em_nuc_sys", data_Em_80_set1);  
 
  
  data_file_580_set1->cd();
  //data_file_580_set1->GetObject("H_Em_nuc_sys", data_Em_580_set1);
  data_file_580_set1->GetObject("H_charge", data_eztar_580_set1); 

  data_file_580_set2->cd();
  //data_file_580_set2->GetObject("H_Em_nuc_sys", data_Em_580_set2);  
  data_file_580_set2->GetObject("H_charge", data_eztar_580_set2); 

  /*
  data_file_750_set1->cd();
  data_file_750_set1->GetObject("H_Em_nuc_sys", data_Em_750_set1);  

  data_file_750_set2->cd();
  data_file_750_set2->GetObject("H_Em_nuc_sys", data_Em_750_set2);  
 
  data_file_750_set3->cd();
  data_file_750_set3->GetObject("H_Em_nuc_sys", data_Em_750_set3);  
  */

  //============================================================================================

  
  TCanvas *c1 = new TCanvas("c1", "", 500, 500);
  data_eztar_580_set1->Add(data_eztar_580_set2);
  data_eztar_580_set1->Draw();
  Double_t err = 0.0;
  Double_t I =  data_eztar_580_set1->IntegralAndError( data_eztar_580_set1->FindBin(-1), data_eztar_580_set1->FindBin(1), err);
  cout << I << " +/- " << err << endl;
  //0.00594638 +/- 0.00198479 counts / mC  (contribution from shms_ztar when -80 < Em < -40 MeV)
  //0.00594638 +/- 0.00198479
  

  
  /*
  TCanvas *c1 = new TCanvas("c1", "", 500, 500);
  data_Em_80_set1->Draw();

  Double_t err_80_dummy = 0.0;
  Double_t err_80_data = 0.0;

  Double_t I_80_dummy = data_Em_80_set1->IntegralAndError(data_Em_80_set1->FindBin(-0.08), data_Em_80_set1->FindBin(-0.02), err_80_dummy);
  Double_t I_80_data = data_Em_80_set1->IntegralAndError(data_Em_80_set1->FindBin(-0.02), data_Em_80_set1->FindBin(0.04), err_80_data);

  cout << "Al_dummy yield = " << I_80_dummy << " +/- " << err_80_dummy << endl;
  cout << "deut+dummy yield = " << I_80_data << " +/- " << err_80_data << endl;
  cout << "fractional dummy yield / (data+dummy) = " <<  (I_80_dummy /  I_80_data) * 100 << endl; 
  
  
  TCanvas *c1 = new TCanvas("c1", "", 500, 500);
  data_Em_580_set1->Add(data_Em_580_set2);
  data_Em_580_set1->Draw();

  Double_t err_580_dummy = 0.0;
  Double_t err_580_data = 0.0;

  Double_t I_580_dummy = data_Em_580_set1->IntegralAndError(data_Em_580_set1->FindBin(-0.08), data_Em_580_set1->FindBin(-0.04), err_580_dummy);
  //Double_t I_580_data = data_Em_580_set1->IntegralAndError(data_Em_580_set1->FindBin(-0.02), data_Em_580_set1->FindBin(0.02), err_580_data);

  cout << "Al_dummy yield = " << I_580_dummy << " +/- " << err_580_dummy << endl;
  //cout << "deut+dummy yield = " << I_580_data << " +/- " << err_580_data << endl;
  //cout << "fractional dummy yield / (data+dummy) = " <<  (I_580_dummy /  I_580_data) * 100 << endl; 
  //Al_dummy yield = 0.00850841 +/- 0.0023623
  //deut+dummy yield = 0.258108 +/- 0.0132694
  //fractional dummy yield / (data+dummy) = 3.29645
  //-------------------------------
  

  
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
  */
}
