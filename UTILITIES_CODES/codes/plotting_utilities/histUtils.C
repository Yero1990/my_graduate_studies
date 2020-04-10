//Utility function that takes as argument 2 histogram objects, and the xlabel, ylabel (assume same range) and plots their ratio
//Used to plot data/simc comparison and their ratio
void hist_ratio(TH1F *hdata, TH1F *hsimc, TString xlabel="", TString ylabel="", TString title="")
{

  gStyle->SetOptStat(0);
  gStyle->SetTitleFontSize(0.08);
  gStyle->SetTitleFont(22, "");
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFont(22);
  gStyle->SetLegendTextSize(0.06);

  //Create Canvas
  auto c = new TCanvas("c", "c", 700, 700);
  c->Divide(1,2, 0., 0.0);

  //VARIABLES TO Normalize histogram (if desired)
  Double_t scale; //Used to scale SIMC histograms by 1./h->Integral(data)
  int dataI;  //data integral
  int simcI;  //simc integral
  
  //Set Histo Aesthetics
  hsimc->SetLineWidth(2);
  hsimc->SetLineColor(kRed);

  hdata->SetFillColorAlpha(kBlue, 0.35);
  hdata->SetFillStyle(3004);

  //Create new histo to store ratio
  TH1F *hratio = (TH1F*) hdata->Clone();

  //Set Histos Axis Label Size
  hdata->GetYaxis()->SetLabelSize(0.06);
  hratio->GetYaxis()->SetLabelSize(0.05);
  hratio->GetXaxis()->SetLabelSize(0.06);
  hdata->SetTitleSize(0.06, "XY");
  hratio->SetTitleSize(0.06, "XY");

  hratio->SetLineWidth(2);
  hratio->SetLineColor(kBlack);

  
  c->cd(1);
  //auto leg = new TLegend(0.1,0.8,0.28,0.9); 
  auto leg = new TLegend(0.14,0.85,0.25,0.65); 

  TPad* pad1 = (TPad*)c->GetPad(1);
  pad1->SetBottomMargin(0.02);
  pad1->SetRightMargin(0.02);
  pad1->SetLeftMargin(0.1);
  pad1->SetTopMargin(0.1);
  pad1->SetFrameLineWidth(2);
  hdata->Draw("samehistE0");
  hsimc->Draw("samesE0");
  hdata->SetTitle(title);

  hdata->GetXaxis()->SetLabelSize(0);
  hdata->GetYaxis()->SetTitle(ylabel);
  hdata->GetYaxis()->CenterTitle();
  hdata->GetYaxis()->SetRangeUser(0,hdata->GetMaximum()+0.6*hdata->GetMaximum());
  hdata->GetYaxis()->SetTitleOffset(0.6);
  hdata->SetLabelFont(22, "XY");
  hdata->SetTitleFont(22, "XY");
  
  dataI = hdata->Integral();
  simcI = hsimc->Integral();
  leg->AddEntry(hdata,Form("Data | Integral: %d", dataI),"f");
  leg->AddEntry(hsimc,Form("SIMC | Integral: %d", simcI));
  //leg->AddEntry((TObject*)0, Form("Data Integral: %d", dataI), "");
  //leg->AddEntry((TObject*)0, Form("SIMC Integral: %d", simcI), "");
  leg->Draw();
  
  c->cd(2);
  TPad* pad = (TPad*)c->GetPad(2);
  pad->SetTopMargin(0.01);
  pad->SetBottomMargin(0.5);
  pad->SetRightMargin(0.02);
  pad->SetLeftMargin(0.1);
  pad->SetFrameLineWidth(2);

  hratio->Divide(hdata, hsimc);
  hratio->GetYaxis()->SetRangeUser(0.,2.1);
  hratio->SetTitle("");

  hratio->GetXaxis()->SetTitle(xlabel);
  hratio->GetYaxis()->SetTitle("Y_{data} / Y_{SIMC}");
  hratio->GetXaxis()->CenterTitle();
  hratio->GetYaxis()->CenterTitle();
  hratio->GetYaxis()->SetTitleOffset(0.6);
  hratio->SetLabelFont(22, "XY");
  hratio->SetTitleFont(22, "XY");

  hratio->Draw();
  pad->Update();
  //Draw lines at +/- 10 % of 1.
  TLine* lmin_10 = new TLine(pad->GetUxmin(),0.9,pad->GetUxmax(),0.9);
  TLine* lmax_10 = new TLine(pad->GetUxmin(),1.1,pad->GetUxmax(),1.1);

  lmin_10->SetLineColor(kBlack);
  lmin_10->SetLineStyle(2);
  lmin_10->Draw();
  lmax_10->SetLineColor(kBlack);
  lmax_10->SetLineStyle(2);
  lmax_10->Draw();
  
  TLine* lmin_20 = new TLine(pad->GetUxmin(),0.8,pad->GetUxmax(),0.8);
  TLine* lmax_20 = new TLine(pad->GetUxmin(),1.2,pad->GetUxmax(),1.2);

  lmin_20->SetLineColor(kBlue);
  lmin_20->SetLineStyle(2);
  lmin_20->Draw();
  lmax_20->SetLineColor(kBlue);
  lmax_20->SetLineStyle(2);
  lmax_20->Draw();
  /*
  auto c1 = new TCanvas("c1", "fit residual simple");
  auto h1 = new TH1D("h1", "h1", 50, -5, 5);
  h1->FillRandom("gaus", 2000);
  h1->Fit("gaus");
  h1->GetXaxis()->SetTitle("x");
  h1->GetYaxis()->SetTitle("y");
  c1->Clear();
  auto rp1 = new TRatioPlot(h1);
  std::vector<double> lines = {-3, -2, -1, 0, 1, 2, 3};
  rp1->SetGridlines(lines);
  rp1->Draw();
  rp1->GetLowerRefGraph()->SetMinimum(-4);
  rp1->GetLowerRefGraph()->SetMaximum(4);
  c1->Update();
  */
}

