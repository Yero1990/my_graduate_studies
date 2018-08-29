void checkREF(Int_t run)
{

  //gROOT->SetBatch(kTRUE); 
  
  TString filename =Form("../../../ROOTfiles/hms_replay_production_all_%d_50000.root", run);                              
      
  //read the file and get the tree                     
  TFile *data_file = new TFile(filename, "READ");                                                     
  TTree *T = (TTree*)data_file->Get("T");

  //Define Plane Names
  TString plane_names[12];
  plane_names[0]="1u1";      //group 81-96 is known to be in the out of sync TDC in HMS crate03 
  plane_names[1]="1u2";   
  plane_names[2]="1x1" ;   
  plane_names[3]="1x2";   
  plane_names[4]="1v1";   
  plane_names[5]="1v2";      //group 49-64 is known to be in the out of sync TDC in HMS crate03   
  plane_names[6]="2v2";    
  plane_names[7]="2v1";     //group 81-96 is known to be in the out of sync TDC in HMS crate03 
  plane_names[8]="2x2";   
  plane_names[9]="2x1";   
  plane_names[10]="2u2";    //group 49-64 is known to be in the out of sync TDC in HMS crate03
  plane_names[11]="2u1"; 
  
  //User Input
  string hreftime = "hDCREF1";    //choose which reference time to use (hDCREF1, 2, 3, 4, or 5)
  string pl_nm = "1v2";   //choose plane name to look at
  Int_t group_min = 49;     //minimum wire in chosen group
  Int_t group_max = 64;      //maximum wire in chosen group

  TCanvas *c1 = new TCanvas(Form("Plane %s, REF Time ", pl_nm.c_str()), Form("Plane %s, Wire Group (%d - %d)", pl_nm.c_str(), group_min, group_max),  1500, 1500);
  c1->Divide(1,2);
 
  //Create Histograms
  TH2F *H_norefCorr = new TH2F("norefCorr", Form("Run %d: hDCRef vs. %s REF. Uncorrected tdcTime: Wire group: %d-%d", run, pl_nm.c_str(), group_min, group_max), 500, 1000, 6000, 500, 1600, 2000);
  TH2F *H_refCorr = new TH2F("refCorr", Form("Run %d: hDCRef vs. %s REF. Corrected tdcTime: Wire group: %d-%d", run, pl_nm.c_str(), group_min, group_max), 500, -16000, -10000, 500, 1660, 2000);

  TCut minCut = Form("H.dc.%s.wirenum>=%d", pl_nm.c_str(), group_min);
  TCut maxCut = Form("H.dc.%s.wirenum<=%d", pl_nm.c_str(), group_max);

  //Draw Histograms
  c1->cd(1);
  T->Draw(Form("T.hms.%s_tdcTime:H.dc.%s.rawnorefcorrtdc>>norefCorr", hreftime.c_str(), pl_nm.c_str()), minCut&&maxCut, "colz");
  c1->cd(2);
  T->Draw(Form("T.hms.%s_tdcTime:H.dc.%s.rawtdc>>refCorr", hreftime.c_str(), pl_nm.c_str()), minCut&&maxCut, "colz");


  c1->SaveAs(Form("./plane%s_%d-%d.pdf", pl_nm.c_str(), group_min, group_max));


}

