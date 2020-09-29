#include "../plotting_utilities/histUtils.C"   //load utilities code for histo plotting

//Script to make comparison between SIMC and Commissioning Data from HallC Spring 2018
//The comparisons are all charge normalized and corrected for all inefficiencies

void add_hist()
{
  
  //gROOT->SetBatch(kTRUE);  
  //gStyle->SetOptStat(1001111);
  gStyle->SetOptStat(0000000);
  string thnq_set;
  string root_dir_1;
  string root_dir_2;

  string root_dir_3;
  string root_dir_4;

  //Deep 580 MeV,  Em_cut (-0.1, -0.04)  -- draw background
  root_dir_1 = "../../root_files/bkg_study_deep_580set1_Em_cut2/deep_data_histos_pm580_set1_combined.root";
  root_dir_2 = "../../root_files/bkg_study_deep_580set2_Em_cut2/deep_data_histos_pm580_set2_combined.root";

  //Deep 580 MeV em_cut (-0.02, 0.04)  -- draw data+background
  root_dir_3 = "../../root_files/bkg_study_deep_580set1_Em_cut1/deep_data_histos_pm580_set1_combined.root";
  root_dir_4 = "../../root_files/bkg_study_deep_580set2_Em_cut1/deep_data_histos_pm580_set2_combined.root";


  //TFile *data_file_80_set1 = new TFile(root_dir_0.c_str());
  
  TFile *data_file_580_set1_bkg = new TFile(root_dir_1.c_str());
  TFile *data_file_580_set2_bkg = new TFile(root_dir_2.c_str());

  TFile *data_file_580_set1 = new TFile(root_dir_3.c_str());
  TFile *data_file_580_set2 = new TFile(root_dir_4.c_str());
  
  TH1F *data_Em_580_set1 = 0;
  TH1F *data_Em_580_set2 = 0;

  TH1F *data_eztar_580_set1_bkg = 0;
  TH1F *data_eztar_580_set2_bkg = 0;

  TH1F *data_eztar_580_set1 = 0;
  TH1F *data_eztar_580_set2 = 0;

  
  //Get DATA (background cuts)
  data_file_580_set1_bkg->cd();
  data_file_580_set1_bkg->GetObject("H_Em_nuc_sys", data_Em_580_set1);

  data_file_580_set1_bkg->GetObject("H_ez_tar", data_eztar_580_set1_bkg); 

  data_file_580_set2_bkg->cd();
  data_file_580_set2_bkg->GetObject("H_Em_nuc_sys", data_Em_580_set2);  
  data_file_580_set2_bkg->GetObject("H_ez_tar", data_eztar_580_set2_bkg); 

  //data+bkg
  data_file_580_set1->cd();
  data_file_580_set1->GetObject("H_ez_tar", data_eztar_580_set1); 
  data_file_580_set2->cd();
  data_file_580_set2->GetObject("H_ez_tar", data_eztar_580_set2); 


  

  //============================================================================================

  TCanvas *c1 = new TCanvas("c1", "", 500, 500);
  c1->cd();

  //Draw Em_nuc_sys (without any cuts)
  data_Em_580_set1->Add(data_Em_580_set2);
  data_Em_580_set1->Draw("E0");
  

  //========================================================

  TCanvas *c2 = new TCanvas("c2", "", 500, 500);
  c2->cd();
  data_eztar_580_set1->Add(data_eztar_580_set2);
  data_eztar_580_set1->Draw("E0");

  data_eztar_580_set1_bkg->Add(data_eztar_580_set2_bkg);
  data_eztar_580_set1_bkg->SetLineColor(kRed);
  data_eztar_580_set1_bkg->Draw("samesE0");
  
  Double_t I_580_data_eztar = 0.0;
  Double_t err_580_data_eztar = 0.0;

  Double_t I_580_bkg_eztar = 0.0;
  Double_t err_580_bkg_eztar = 0.0;
  
  I_580_data_eztar = data_eztar_580_set1->IntegralAndError(data_eztar_580_set1->FindBin(-7), data_eztar_580_set1->FindBin(7), err_580_data_eztar );
  I_580_bkg_eztar = data_eztar_580_set1_bkg->IntegralAndError(data_eztar_580_set1_bkg->FindBin(-7), data_eztar_580_set1_bkg->FindBin(7), err_580_bkg_eztar );
  
  cout << "data+bkg eztar yield = " <<  I_580_data_eztar << " +/- " << err_580_data_eztar << endl;
  cout << "bkg eztar yield = " <<  I_580_bkg_eztar << " +/- " << err_580_bkg_eztar << endl;
  cout << "bkg / (data+bkg) = " <<  (I_580_bkg_eztar / I_580_data_eztar) * 100. << " % " << endl;

  
}
