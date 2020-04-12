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
  TCanvas *c = new TCanvas(title, "c", 700, 700);
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
  hdata->SetLineWidth(2);
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
  auto leg = new TLegend(0.14,0.84,0.25,0.64); 

  TPad* pad1 = (TPad*)c->GetPad(1);
  pad1->SetBottomMargin(0.02);
  pad1->SetRightMargin(0.02);
  pad1->SetLeftMargin(0.1);
  pad1->SetTopMargin(0.13);
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

  double dataI_err, simcI_err;
  double nbins = hdata->GetNbinsX();  //Get total number of bins (excluding overflow)
  dataI = hdata->IntegralAndError(1, nbins, dataI_err);
  simcI = hsimc->IntegralAndError(1, nbins, simcI_err);
  double R = (float)dataI / simcI;
  double R_err = R * sqrt(pow(dataI_err/dataI,2) + pow(simcI_err/simcI,2));
  
  leg->AddEntry(hdata,Form("Data | Integral: %d", dataI),"f");
  leg->AddEntry(hsimc,Form("SIMC | Integral: %d", simcI));
  leg->AddEntry((TObject*)0, Form("Ratio: %.3f #pm %.3f", R, R_err), "");
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
}

void plot_hist(TH1F *hist, TString xlabel="", TString ylabel="", TString title="", TString set_logy="")
{

  gStyle->SetOptStat(0);
  gStyle->SetTitleFontSize(0.05);
  gStyle->SetTitleFont(22, "");
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFont(22);
  gStyle->SetLegendTextSize(0.03);

  //Create Canvas
  TCanvas *c = new TCanvas("c", "c", 900, 700);
  c->SetFrameLineWidth(2);
  //VARIABLES TO Normalize histogram (if desired)
  Double_t scale; //Used to scale SIMC histograms by 1./h->Integral(data)
  int dataI;  //data, simc integral
  double dataI_err;


  hist->SetFillColorAlpha(kBlue, 0.35);
  hist->SetFillStyle(3004);

  //Set Histos Axis Label Size
  hist->GetYaxis()->SetLabelSize(0.06);
  hist->SetTitleSize(0.04, "XY");

  TLegend *leg = new TLegend(0.14,0.8,0.4,0.7); 

  if(set_logy=="logy"){
    c->SetLogy();
  }

  hist->Draw("samehistE0");

  
  hist->SetTitle(title);
  
  hist->GetXaxis()->SetLabelSize(0.04);
  hist->GetYaxis()->SetLabelSize(0.04);
  
  hist->GetYaxis()->SetTitle(ylabel);
  hist->GetXaxis()->SetTitle(xlabel);

  hist->GetYaxis()->CenterTitle();
  hist->GetXaxis()->CenterTitle();

  
  //if(set_logy==""){hist->GetYaxis()->SetRangeUser(0,hist->GetMaximum()+0.6*hist->GetMaximum());}

  hist->GetYaxis()->SetTitleOffset(1.3);
  hist->GetXaxis()->SetTitleOffset(1.3);
  
  hist->SetLabelFont(22, "XY");
  hist->SetTitleFont(22, "XY");

  double nbins = hist->GetNbinsX();  //Get total number of bins (excluding overflow)
  dataI = hist->IntegralAndError(1, nbins, dataI_err);
  
  leg->AddEntry(hist,Form("Data | Integral: %d", dataI),"f");
  leg->Draw();

}



//---------------------

void compare_hist(TH1F *hdata, TH1F *hsimc, TString xlabel="", TString ylabel="", TString title="")
{

  gStyle->SetOptStat(0);
  gStyle->SetTitleFontSize(0.04);
  gStyle->SetTitleFont(22, "");
  gStyle->SetLegendBorderSize(0);
  gStyle->SetLegendFont(22);
  gStyle->SetLegendTextSize(0.03);

  //Create Canvas
  TCanvas *c = new TCanvas(title, "c", 900, 700);
  c->SetFrameLineWidth(2);

  //VARIABLES TO Normalize histogram (if desired)
  Double_t scale; //Used to scale SIMC histograms by 1./h->Integral(data)
  int dataI;  //data integral
  int simcI;  //simc integral
  
  //Set Histo Aesthetics
  hsimc->SetLineWidth(2);
  hsimc->SetLineColor(kRed);

  hdata->SetFillColorAlpha(kBlue, 0.35);
  hdata->SetFillStyle(3004);

  //Set Histos Axis Label Size
  hdata->GetYaxis()->SetLabelSize(0.04);
  hdata->GetXaxis()->SetLabelSize(0.04);
  hdata->SetTitleSize(0.04, "XY");

 
  //auto leg = new TLegend(0.1,0.8,0.28,0.9); 
  TLegend *leg = new TLegend(0.14,0.88,0.25,0.73); 

  hdata->Draw("samehistE0");
  hsimc->Draw("samesE0");
  hdata->SetTitle(title);

  hdata->GetYaxis()->SetTitle(ylabel);
  hdata->GetXaxis()->SetTitle(xlabel);
  hdata->GetYaxis()->CenterTitle();
  hdata->GetXaxis()->CenterTitle();
  hdata->GetYaxis()->SetRangeUser(0,hdata->GetMaximum()+0.6*hdata->GetMaximum());
  hdata->GetXaxis()->SetTitleOffset(1.);  
  hdata->GetYaxis()->SetTitleOffset(1.);
  hdata->SetLabelFont(22, "XY");
  hdata->SetTitleFont(22, "XY");

  double dataI_err, simcI_err;
  double nbins = hdata->GetNbinsX();  //Get total number of bins (excluding overflow)
  dataI = hdata->IntegralAndError(1, nbins, dataI_err);
  simcI = hsimc->IntegralAndError(1, nbins, simcI_err);
  double R = (float)dataI / simcI;
  double R_err = R * sqrt(pow(dataI_err/dataI,2) + pow(simcI_err/simcI,2));
  
  leg->AddEntry(hdata,Form("Data | Integral: %d", dataI),"f");
  leg->AddEntry(hsimc,Form("SIMC | Integral: %d", simcI));
  leg->AddEntry((TObject*)0, Form("Ratio: %.3f #pm %.3f", R, R_err), "");
  leg->Draw();


}
