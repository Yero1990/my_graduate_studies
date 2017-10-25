void plot_eff()
{
//Plot Efficiency Curve for HV scan of DC2, at a threshold = 4.5 V:
  //HMS Runs 673 (HV = -1940, 435 evts) 
  //         674 (HV = -1840, 557 evts)
  //         675 (HV = -1740, 554 evts)
  //         676 (HV = -1640, 523 evts)   //not enough statistics
  //         677 (HV = -1970, 547 evts)
  //         678 (HV = -1700, 1294 evts)
  //         683 (HV = -1900, 725 evts )
  // create first graph

  const Int_t n = 6;
  
  //Read in the data
  Double_t HV[] = {1700, 1740, 1840, 1900, 1940, 1970};   //X-axis (High Voltage)
  Double_t ex[] = {0.0,0.0,0.0,0.0,0.0,0.0};            // no error in X
  
  //HMS DC2 Plane Efficiency
  Double_t eff_u1[] = {0.72, 0.79  , 0.974, 0.998, 1.0, 1.0};               
  Double_t err_u1[] = {0.06, 0.0452, 0.0069,0.00139, 0.0, 0.0};    //uncertainty is given in percent    
  
  Double_t eff_u2[] = {0.72, 0.853 , 0.986, 0.998,  1.0, 1.0};            
  Double_t err_u2[] = {0.06, 0.0408, 0.0052, 0.00139, 0.0, 0.0};      

  Double_t eff_x1[] = {0.8, 0.771 , 0.978, 0.998, 1.0, 1.0};               
  Double_t err_x1[] = {0.059, 0.0461, 0.0064, 0.00139, 0.0, 0.0};      
  
  Double_t eff_x2[] = {0.8, 0.888,  0.968, 0.9958, 0.995, 0.994};            
  Double_t err_x2[] = {0.059, 0.037,  0.00768, 0.00241, 0.00328, 0.00319};      

  Double_t eff_v1[] = {0.923, 0.901,  0.990, 1.0,  1.0, 1.0};               
  Double_t err_v1[] = {0.042, 0.0353, 0.0044,0.0, 0.0, 0.0};      
  
  Double_t eff_v2[] = {0.9, 0.864,  1.0, 1.0, 1.0, 1.0};            
  Double_t err_v2[] = {0.047, 0.0397, 0.0, 0.0, 0.0, 0.0};      


  TCanvas *c1 = new TCanvas("c1", "", 2000,500);
  c1->SetGrid();

  // draw a frame to define the range
  TMultiGraph *mg = new TMultiGraph();
  
  TGraphErrors *gr1 = new TGraphErrors(n,HV,eff_u1,ex,err_u1);
  gr1->SetMarkerColor(kBlue);
  gr1->SetMarkerStyle(20);
  mg->Add(gr1);
  
  TGraphErrors *gr2 = new TGraphErrors(n,HV,eff_u2,ex,err_u2);
  gr2->SetMarkerColor(kRed);
  gr2->SetMarkerStyle(20);
  mg->Add(gr2);

  TGraphErrors *gr3 = new TGraphErrors(n,HV,eff_x1,ex,err_x1);
  gr3->SetMarkerColor(kBlack);
  gr3->SetMarkerStyle(21);
  mg->Add(gr3);
  
  TGraphErrors *gr4 = new TGraphErrors(n,HV,eff_x2,ex,err_x2);
  gr4->SetMarkerColor(kMagenta);
  gr4->SetMarkerStyle(21);
  mg->Add(gr4);

   TGraphErrors *gr5 = new TGraphErrors(n,HV,eff_v1,ex,err_v1);
  gr5->SetMarkerColor(8);
  gr5->SetMarkerStyle(22);
  mg->Add(gr5);
  
  TGraphErrors *gr6 = new TGraphErrors(n,HV,eff_v2,ex,err_v2);
  gr6->SetMarkerColor(kCyan);
  gr6->SetMarkerStyle(22);
  mg->Add(gr6);

  mg->Draw("ap");
  mg->SetTitle("HMS DCII High Voltage Scan at Threshold = 4.5 V");
  mg->GetXaxis()->SetTitle("High Voltage (Volts)");
  mg->GetYaxis()->SetTitle("Plane Efficiency (%)");
  mg->GetXaxis()->CenterTitle();
  mg->GetYaxis()->CenterTitle();
  
  //Add a legend
  auto legend = new TLegend(0.1, 0.2, 0.48, 0.9);
  legend->SetTextFont(0);
  legend->SetHeader("HMS DCII Planes", "C");
  legend->AddEntry(gr1,"Plane U ", "lep");
  legend->AddEntry(gr2,"Plane U'", "lep");
  legend->AddEntry(gr3,"Plane X", "lep");
  legend->AddEntry(gr4,"Plane X'", "lep");
  legend->AddEntry(gr5,"Plane V", "lep");
  legend->AddEntry(gr6,"Plane V'", "lep");	
  legend->Draw();

}
