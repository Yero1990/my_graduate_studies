// 10/03/2017. here I want to compare the dark current drawn by hms ( with ARCO2) shms(ArCO2) and hms (ArC2H5)



#include <iostream>
#include <string>
#include "TGraph.h"
#include "TCanvas.h"
using namespace std;
void chamber()
{
  gStyle->SetOptStat(0);
  TCanvas *c1=new TCanvas("c1","",1800,1800);
  c1->SetGrid();
  TMultiGraph *mg = new TMultiGraph();
  c1->cd();

  // Hms with ArCO2
  mg->SetTitle("Comparison of Current Drawn with Different Gas Mixtures");
  Double_t x1[7]={1750,1775,1800,1825,1850,1875,1900};
  Double_t y1[7] ={3.9,4.2,4.7,5.5,7.2,29,105};
  TGraph *gr1 = new TGraph(7,x1,y1);
  gr1->SetMarkerStyle(kFullCircle);
  gr1->SetMarkerColor(4);
  mg->Add(gr1); 
  // SHMS ArCO2
  Double_t x2[9]={1750,1800,1825,1850,1875,1900,1925,1950,1962};
  Double_t y2[9] ={0.7,0.9,1.0,1.1,1.3,1.5,2.3,18,45};
  TGraph *gr2 = new TGraph(9,x2,y2);
  gr2->SetMarkerStyle(kFullSquare);
  gr2->SetMarkerColor(kBlack);
  mg->Add(gr2);
  //HMS Ar:Ethane 50:50
  Double_t x3[17]={1800,1825,1850,1875,1900,1925,1950,1975,2000,2025,2050,2075,2100,2125,2150,2175,2200};
  Double_t y3[17] ={0.015,0.02,0.02,0.02,0.02,0.025,0.03,0.035,0.04,0.065,0.085,0.1,0.11,0.125,0.17,0.3,21.1};
 TGraph *gr3 = new TGraph(17,x3,y3);

 gr3->SetMarkerStyle(kFullTriangleUp);
  gr3->SetMarkerColor(2);
  mg->Add(gr3);
  mg->Draw("alp");

  
  auto leg = new TLegend(0.1, 0.7, 0.3, 0.9);
  // leg->SetTextFont(0);
   // leg->SetFillColor(0);
  leg->SetHeader("Gas Mixtures");
  leg->AddEntry(gr1, "HMS DCI, 75:25 Ar/CO_{2}", "lp");
  leg->AddEntry(gr2, "SHMS DCI, 75:25 Ar/CO_{2} ", "lp");
  leg->AddEntry(gr3, "HMS DCII, 50:50 Ar/Ethane", "lp");
  leg->Draw();
  
  int m = 1000;
  //HMS Ar:Ethane 50:50 (Current in nanoamps to put in subplot)
  Double_t x4[17]={1800,1825,1850,1875,1900,1925,1950,1975,2000,2025,2050,2075,2100,2125,2150,2175,2200};
  Double_t y4[17] ={0.015*m,0.02*m,0.02*m,0.02*m,0.02*m,0.025*m,0.03*m,0.035*m,0.04*m,0.065*m,0.085*m,0.1*m,0.11*m,0.125*m,0.17*m,0.3*m,21.1*m};
 TGraph *gr4 = new TGraph(17,x4,y4);

 gr4->SetMarkerStyle(kFullTriangleUp);
  gr4->SetMarkerColor(2);

  TPad *subpad = new TPad("subpad","",0.5,0.5,0.88,0.88);
  subpad->Draw();
  subpad->cd();
  gr4->SetTitle("HMS DCII");
  gr4->Draw();
  

  
}
