void histproj()
{
  TString fname1 = "deep_data_histos_pm80_set1_combined.root";
  TString fname2 = "deep_data_histos_pm80_set1_combined_radcorr.root";
  TString fname3 = "deep_simc_histos_pm80_lagetfsi_RadCorrRatio_set1.root";

  TFile *file1 = new TFile(fname1);
  TFile *file2 = new TFile(fname2);
  TFile *file3 = new TFile(fname3);

  
  //Before rad. corr
  file1->cd();
  TH2F *h1 = 0;
  file1->GetObject("H_Pm_vs_thnq", h1);
  //h1->Draw("colz");
  TH1D *histo1= h1->ProjectionY("name1",h1->GetXaxis()->FindBin(30.1),h1->GetXaxis()->FindBin(39.9),"");
  histo1->Draw();

  
  //After rad corr
  file2->cd();
  TH2F *h2 = 0;
  file2->GetObject("H_Pm_vs_thnq", h2);
  //h2->Draw("colz");
  TH1D *histo2= h2->ProjectionY("name2",h2->GetXaxis()->FindBin(30.1),h2->GetXaxis()->FindBin(39.9),"");
  histo2->Draw("sames");
  
  /*
  //rad corr. ratio
  file3->cd();
  TH2F *h3 = 0;
  file3->GetObject("H_Pm_vs_thnq_ratio", h3);
  //h3->Draw("colz");
  TH1D *histo3= h3->ProjectionY("name2",h3->GetXaxis()->FindBin(30.1),h3->GetXaxis()->FindBin(39.9),"");
  histo3->Draw();
  */
}
