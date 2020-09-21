//Script to make comparison between SIMC and Commissioning Data from HallC Spring 2018

// The SIMC histograms are normalized to the area of the data histograms
// The ratio of data/SIMC is also taken.
// This is to check spectrometer acceptances

//Compare Target Reconstruction/FOCAL PLANE/ Kinematics Variables

void ratioPlot(int run)
{

  int ratio = 1;
  //gROOT->SetBatch(kTRUE);  
  //gStyle->SetOptStat(1001111);
  gStyle->SetOptStat(0000000);


  TString root_dir = "../../root_files/pre_HEEP_ELASTICS/";
  //Pre-defined SIMC/data root file names containing histogram object to comapare
  TString simc_filename =  Form("heep_simc_histos_%d_rad.root", run);                      
  TString data_filename = Form("heep_data_histos_%d_combined.root",run); 

  TString simc_f = root_dir + simc_filename;
  TString data_f = root_dir + data_filename;
  
  //Where to store plots
  string plots_dir = "./";
  string plots_path;
  
  //Open SIMC/data ROOT files;
  TFile *simc_file = new TFile(simc_f);
  TFile *data_file = new TFile(data_f);


  
  //---------------Target ----------------
  //Define SIMC histos ('h'-->hadron arm,  'e'-->electron arm)
  
  TH1F *simc_xtarH =  0;
  TH1F *simc_ytarH =  0;
  TH1F *simc_ztarH =  0;

  TH1F *simc_xtarP = 0;
  TH1F *simc_ytarP =  0;
  TH1F *simc_ztarP =  0;  

  //Define data histos
  TH1F *data_xtarH = 0;
  TH1F *data_ytarH = 0;
  TH1F *data_ztarH = 0;

  TH1F *data_xtarP = 0;                                                                                                     
  TH1F *data_ytarP = 0;                                                                                                                                   
  TH1F *data_ztarP = 0; 

  //---------------Target Reconstruction Variables----------------
  //Define SIMC histos ('h'-->hadron arm,  'e'-->electron arm)
  TH1F *simc_eytar =  0;
  TH1F *simc_exptar =  0;
  TH1F *simc_eyptar =  0;
  TH1F *simc_edelta =  0;

  TH1F *simc_hytar =  0;
  TH1F *simc_hxptar =  0;
  TH1F *simc_hyptar =  0;
  TH1F *simc_hdelta =  0;

  //Define data histos
  TH1F *data_eytar = 0;
  TH1F *data_exptar =  0;
  TH1F *data_eyptar =  0;
  TH1F *data_edelta =  0;

  TH1F *data_hytar = 0;
  TH1F *data_hxptar =  0;
  TH1F *data_hyptar =  0;
  TH1F *data_hdelta =  0;

  //-----------------------------------------------------------
 
  //--------------FOCAL PLANE VARIABLES------------------------

 //Define SIMC histos ('h'-->hadron arm,  'e'-->electron arm)
  TH1F *simc_exfp =  0;
  TH1F *simc_eyfp =  0;
  TH1F *simc_expfp =  0;
  TH1F *simc_eypfp =  0;

  TH1F *simc_hxfp =  0;
  TH1F *simc_hyfp =  0;
  TH1F *simc_hxpfp =  0;
  TH1F *simc_hypfp =  0;
  
  //Define data histos
  TH1F *data_exfp =  0;
  TH1F *data_eyfp =  0;
  TH1F *data_expfp =  0;
  TH1F *data_eypfp =  0;

  TH1F *data_hxfp =  0;
  TH1F *data_hyfp =  0;
  TH1F *data_hxpfp =  0;
  TH1F *data_hypfp =  0;

  //--------------------------------------------------------------

  //-------------------------KINEMATICS---------------------------
  TH1F *simc_Q2 =  0;
  TH1F *simc_omega =  0;
  TH1F *simc_W =  0;
  TH1F *simc_thq = 0;

  TH1F *simc_xbj = 0;
  TH1F *simc_th_elec = 0;                                  
  TH1F *simc_kf = 0;  
  TH1F *simc_emiss = 0;

  //Kinematics 2
  TH1F *simc_Pm = 0;
  TH1F *simc_Pf = 0;
  TH1F *simc_th_prot = 0;
  TH1F *simc_q = 0;    //q-vector magnitude
  TH1F *simc_thpq = 0;
  TH1F *simc_Pmx = 0;
  TH1F *simc_Pmy = 0;
  TH1F *simc_Pmz = 0;

  //Define data histos
  TH1F *data_Q2 =  0;
  TH1F *data_omega =  0;
  TH1F *data_W =  0;
  TH1F *data_thq = 0;

  TH1F *data_xbj = 0;
  TH1F *data_th_elec = 0;
  TH1F *data_kf = 0;
  TH1F *data_emiss = 0;

   //Kinematics 2
  TH1F *data_Pm = 0;
  TH1F *data_Pf = 0;
  TH1F *data_th_prot = 0;
  TH1F *data_q = 0;    //q-vector magnitude
  TH1F *data_thpq = 0;
  TH1F *data_Pmx = 0;
  TH1F *data_Pmy = 0;
  TH1F *data_Pmz = 0;

  //---------------------------------------------------------------

 //change to simc_file
  simc_file->cd();

  
  //----------Get Target Histograms------------------
  //Get Histogram objects from SIMC rootfile
  simc_file->GetObject("H_hx_tar", simc_xtarH);
  simc_file->GetObject("H_hy_tar", simc_ytarH);
  simc_file->GetObject("H_hz_tar", simc_ztarH);

  simc_file->GetObject("H_ex_tar", simc_xtarP);
  simc_file->GetObject("H_ey_tar", simc_ytarP);  
  simc_file->GetObject("H_ez_tar", simc_ztarP); 

  //Set SIMC Histo Aesthetics
  simc_xtarH->SetLineColor(kRed);
  simc_xtarH->SetLineWidth(2);
  simc_ytarH->SetLineColor(kRed);
  simc_ytarH->SetLineWidth(2);
  simc_ztarH->SetLineColor(kRed);
  simc_ztarH->SetLineWidth(2);

  simc_xtarP->SetLineColor(kRed);
  simc_xtarP->SetLineWidth(2);
  simc_ytarP->SetLineColor(kRed);          
  simc_ytarP->SetLineWidth(2);           
  simc_ztarP->SetLineColor(kRed);                        
  simc_ztarP->SetLineWidth(2);  
  
  //change to data_file
  data_file->cd();

  //Get Histogram objects from data rootfile
  data_file->GetObject("H_hx_tar", data_xtarH);
  data_file->GetObject("H_hy_tar", data_ytarH);
  data_file->GetObject("H_hz_tar", data_ztarH);

  data_file->GetObject("H_ex_tar", data_xtarP); 
  data_file->GetObject("H_ey_tar", data_ytarP);                        
  data_file->GetObject("H_ez_tar", data_ztarP); 
  
  //Set data Histo Aesthetics
  data_xtarH->SetFillColorAlpha(kBlue, 0.35);
  data_xtarH->SetFillStyle(3004);
  data_ytarH->SetFillColorAlpha(kBlue, 0.35);
  data_ytarH->SetFillStyle(3004);
  data_ztarH->SetFillColorAlpha(kBlue, 0.35);
  data_ztarH->SetFillStyle(3004);

  data_xtarP->SetFillColorAlpha(kBlue, 0.35);         
  data_xtarP->SetFillStyle(3004); 
  data_ytarP->SetFillColorAlpha(kBlue, 0.35);                                  
  data_ytarP->SetFillStyle(3004);                               
  data_ztarP->SetFillColorAlpha(kBlue, 0.35);             
  data_ztarP->SetFillStyle(3004);  

  //-----------------------------------------------------------------


  //---------------------------------------------------------------

 //change to simc_file
  simc_file->cd();

  //----------Get Target Reconstructed Histograms------------------
  //Get Histogram objects from SIMC rootfile
  simc_file->GetObject("H_eytar", simc_eytar);
  simc_file->GetObject("H_exptar", simc_exptar);
  simc_file->GetObject("H_eyptar", simc_eyptar);
  simc_file->GetObject("H_edelta", simc_edelta);

  simc_file->GetObject("H_hytar", simc_hytar);
  simc_file->GetObject("H_hxptar", simc_hxptar);
  simc_file->GetObject("H_hyptar", simc_hyptar);
  simc_file->GetObject("H_hdelta", simc_hdelta);

  simc_file->GetObject("H_Pmx_Lab", simc_Pmx);
  simc_file->GetObject("H_Pmy_Lab", simc_Pmy);
  simc_file->GetObject("H_Pmz_Lab", simc_Pmz);

  //Set SIMC Histo Aesthetics
  simc_eytar->SetLineColor(kRed);
  simc_eytar->SetLineWidth(2);
  simc_exptar->SetLineColor(kRed);
  simc_exptar->SetLineWidth(2);
  simc_eyptar->SetLineColor(kRed);
  simc_eyptar->SetLineWidth(2);
  simc_edelta->SetLineColor(kRed);
  simc_edelta->SetLineWidth(2);
  
  simc_hytar->SetLineColor(kRed);
  simc_hytar->SetLineWidth(2);
  simc_hxptar->SetLineColor(kRed);
  simc_hxptar->SetLineWidth(2);
  simc_hyptar->SetLineColor(kRed);
  simc_hyptar->SetLineWidth(2);
  simc_hdelta->SetLineColor(kRed);
  simc_hdelta->SetLineWidth(2);

  simc_Pmx->SetLineColor(kRed);
  simc_Pmx->SetLineWidth(2);
  simc_Pmy->SetLineColor(kRed);
  simc_Pmy->SetLineWidth(2);
  simc_Pmz->SetLineColor(kRed);
  simc_Pmz->SetLineWidth(2);

  //change to data_file
  data_file->cd();

  //Get Histogram objects from data rootfile
  data_file->GetObject("H_eytar", data_eytar);
  data_file->GetObject("H_exptar", data_exptar);
  data_file->GetObject("H_eyptar", data_eyptar);
  data_file->GetObject("H_edelta", data_edelta);
  
  data_file->GetObject("H_hytar", data_hytar);
  data_file->GetObject("H_hxptar", data_hxptar);
  data_file->GetObject("H_hyptar", data_hyptar);
  data_file->GetObject("H_hdelta", data_hdelta);

  //Set data Histo Aesthetics
  data_eytar->SetFillColorAlpha(kBlue, 0.35);
  data_eytar->SetFillStyle(3004);
  data_exptar->SetFillColorAlpha(kBlue, 0.35);
  data_exptar->SetFillStyle(3004);
  data_eyptar->SetFillColorAlpha(kBlue, 0.35);
  data_eyptar->SetFillStyle(3004);
  data_edelta->SetFillColorAlpha(kBlue, 0.35);
  data_edelta->SetFillStyle(3004);

  data_hytar->SetFillColorAlpha(kBlue, 0.35);
  data_hytar->SetFillStyle(3004);
  data_hxptar->SetFillColorAlpha(kBlue, 0.35);
  data_hxptar->SetFillStyle(3004);
  data_hyptar->SetFillColorAlpha(kBlue, 0.35);
  data_hyptar->SetFillStyle(3004);
  data_hdelta->SetFillColorAlpha(kBlue, 0.35);
  data_hdelta->SetFillStyle(3004);

  //-----------------------------------------------------------------

  

  //---------------Get FOCAL PLANE Histograms------------------------

   //change to simc_file
  simc_file->cd();

  //Get Histogram objects from SIMC rootfile
  simc_file->GetObject("H_exfp", simc_exfp);
  simc_file->GetObject("H_eyfp", simc_eyfp);
  simc_file->GetObject("H_expfp", simc_expfp);
  simc_file->GetObject("H_eypfp", simc_eypfp);

  simc_file->GetObject("H_hxfp", simc_hxfp);
  simc_file->GetObject("H_hyfp", simc_hyfp);
  simc_file->GetObject("H_hxpfp", simc_hxpfp);
  simc_file->GetObject("H_hypfp", simc_hypfp);
  //Set SIMC Histo Aesthetics
  simc_exfp->SetLineColor(kRed);
  simc_exfp->SetLineWidth(2);
  simc_eyfp->SetLineColor(kRed);
  simc_eyfp->SetLineWidth(2);
  simc_expfp->SetLineColor(kRed);
  simc_expfp->SetLineWidth(2);
  simc_eypfp->SetLineColor(kRed);
  simc_eypfp->SetLineWidth(2);
  
  simc_hxfp->SetLineColor(kRed);
  simc_hxfp->SetLineWidth(2);
  simc_hyfp->SetLineColor(kRed);
  simc_hyfp->SetLineWidth(2);
  simc_hxpfp->SetLineColor(kRed);
  simc_hxpfp->SetLineWidth(2);
  simc_hypfp->SetLineColor(kRed);
  simc_hypfp->SetLineWidth(2);
  //change to data_file
  data_file->cd();

  //Get Histogram objects from data rootfile
  data_file->GetObject("H_exfp", data_exfp);
  data_file->GetObject("H_eyfp", data_eyfp);
  data_file->GetObject("H_expfp", data_expfp);
  data_file->GetObject("H_eypfp", data_eypfp);

  data_file->GetObject("H_hxfp", data_hxfp);
  data_file->GetObject("H_hyfp", data_hyfp);
  data_file->GetObject("H_hxpfp", data_hxpfp);
  data_file->GetObject("H_hypfp", data_hypfp);
  //Set data Histo Aesthetics
  data_exfp->SetFillColorAlpha(kBlue, 0.35);
  data_exfp->SetFillStyle(3004);
  data_eyfp->SetFillColorAlpha(kBlue, 0.35);
  data_eyfp->SetFillStyle(3004);
  data_expfp->SetFillColorAlpha(kBlue, 0.35);
  data_expfp->SetFillStyle(3004);
  data_eypfp->SetFillColorAlpha(kBlue, 0.35);
  data_eypfp->SetFillStyle(3004);

  data_hxfp->SetFillColorAlpha(kBlue, 0.35);
  data_hxfp->SetFillStyle(3004);
  data_hyfp->SetFillColorAlpha(kBlue, 0.35);
  data_hyfp->SetFillStyle(3004);
  data_hxpfp->SetFillColorAlpha(kBlue, 0.35);
  data_hxpfp->SetFillStyle(3004);
  data_hypfp->SetFillColorAlpha(kBlue, 0.35);
  data_hypfp->SetFillStyle(3004);

  //--------------------------------------------------------------
  
  //------------------Get KINEMATICS VARIABLES--------------------

   //change to simc_file
  simc_file->cd();

  //Get Histogram objects from SIMC rootfile
  simc_file->GetObject("H_Q2", simc_Q2);
  simc_file->GetObject("H_omega", simc_omega);
  simc_file->GetObject("H_W", simc_W);
  simc_file->GetObject("H_theta_q", simc_thq);

  simc_file->GetObject("H_xbj", simc_xbj);
  simc_file->GetObject("H_theta_elec", simc_th_elec);
  simc_file->GetObject("H_kf", simc_kf);
  simc_file->GetObject("H_Emiss", simc_emiss);

  simc_file->GetObject("H_Pm", simc_Pm);
  simc_file->GetObject("H_Pf", simc_Pf);
  simc_file->GetObject("H_theta_prot", simc_th_prot);
  simc_file->GetObject("H_q", simc_q);
  simc_file->GetObject("H_theta_pq", simc_thpq);

  
  
  //Set SIMC Histo Aesthetics
  simc_Q2->SetLineColor(kRed);
  simc_Q2->SetLineWidth(2);
  simc_omega->SetLineColor(kRed);
  simc_omega->SetLineWidth(2);
  simc_W->SetLineColor(kRed);
  simc_W->SetLineWidth(2);
  simc_thq->SetLineColor(kRed);
  simc_thq->SetLineWidth(2);
  
  simc_xbj->SetLineColor(kRed);
  simc_xbj->SetLineWidth(2);
  simc_th_elec->SetLineColor(kRed);
  simc_th_elec->SetLineWidth(2);
  simc_kf->SetLineColor(kRed);
  simc_kf->SetLineWidth(2);
  simc_emiss->SetLineColor(kRed);
  simc_emiss->SetLineWidth(2);
  
  simc_Pm->SetLineColor(kRed);
  simc_Pm->SetLineWidth(2);
  simc_Pf->SetLineColor(kRed);
  simc_Pf->SetLineWidth(2);
  simc_th_prot->SetLineColor(kRed);
  simc_th_prot->SetLineWidth(2);
  simc_q->SetLineColor(kRed);
  simc_q->SetLineWidth(2);
  simc_thpq->SetLineColor(kRed);
  simc_thpq->SetLineWidth(2);

  //change to data_file
  data_file->cd();
  
  //Get Histogram objects from data rootfile
  data_file->GetObject("H_Q2", data_Q2);
  data_file->GetObject("H_omega", data_omega);
  data_file->GetObject("H_W", data_W);
  data_file->GetObject("H_theta_q", data_thq);

  
  data_file->GetObject("H_xbj", data_xbj);
  data_file->GetObject("H_theta_elec", data_th_elec);
  data_file->GetObject("H_kf", data_kf);
  data_file->GetObject("H_Emiss", data_emiss);

  data_file->GetObject("H_Pm", data_Pm);
  data_file->GetObject("H_Pf", data_Pf);
  data_file->GetObject("H_theta_prot", data_th_prot);
  data_file->GetObject("H_q", data_q);
  data_file->GetObject("H_theta_pq", data_thpq);

  data_file->GetObject("H_Pmx_Lab", data_Pmx);
  data_file->GetObject("H_Pmy_Lab", data_Pmy);
  data_file->GetObject("H_Pmz_Lab", data_Pmz);

  //Set data Histo Aesthetics
  data_Q2->SetFillColorAlpha(kBlue, 0.35);
  data_Q2->SetFillStyle(3004);
  data_omega->SetFillColorAlpha(kBlue, 0.35);
  data_omega->SetFillStyle(3004);
  data_W->SetFillColorAlpha(kBlue, 0.35);
  data_W->SetFillStyle(3004);
  data_thq->SetFillColorAlpha(kBlue, 0.35);
  data_thq->SetFillStyle(3004);

  data_xbj->SetFillColorAlpha(kBlue, 0.35);
  data_xbj->SetFillStyle(3004);
  data_th_elec->SetFillColorAlpha(kBlue, 0.35);
  data_th_elec->SetFillStyle(3004);
  data_kf->SetFillColorAlpha(kBlue, 0.35);
  data_kf->SetFillStyle(3004);
  data_emiss->SetFillColorAlpha(kBlue, 0.35);
  data_emiss->SetFillStyle(3004);

  data_Pm->SetFillColorAlpha(kBlue,0.35);
  data_Pm->SetFillStyle(3004);
  data_Pf->SetFillColorAlpha(kBlue,0.35);
  data_Pf->SetFillStyle(3004);
  data_th_prot->SetFillColorAlpha(kBlue,0.35);
  data_th_prot->SetFillStyle(3004);
  data_q->SetFillColorAlpha(kBlue,0.35);
  data_q->SetFillStyle(3004);
  data_thpq->SetFillColorAlpha(kBlue,0.35);
  data_thpq->SetFillStyle(3004);

  data_Pmx->SetFillColorAlpha(kBlue,0.35);
  data_Pmx->SetFillStyle(3004);
  data_Pmy->SetFillColorAlpha(kBlue,0.35);
  data_Pmy->SetFillStyle(3004);
  data_Pmz->SetFillColorAlpha(kBlue,0.35);
  data_Pmz->SetFillStyle(3004);



  //NORMALIZATION VARIABLES
  Double_t scale; //Used to scale SIMC histograms by 1./h->Integral(data)
  int dataI;
  int simcI;
  
  //Overlay SIMC/data plots (*** VERY IMPORTANT ***: Range and #bins must be same)
   //Set Legend 
   auto leg5 = new TLegend(0.1,0.8,0.28,0.9); 
   auto leg6 = new TLegend(0.1,0.8,0.28,0.9);
   auto leg7 = new TLegend(0.1,0.8,0.28,0.9);
   auto leg8 = new TLegend(0.1,0.8,0.28,0.9);


   
   //-----------------PLOT Target Reconstructed Variables SIMC/Data comparison-----------------------

   //Create A Canvas to store Target Recon. variable comparisons in HADRON ARM
   
   TCanvas *c1 = new TCanvas("c1", "Electron Arm: Target Reconstruction", 2000, 1000);
   c1->Divide(2,2);

   c1->cd(1);
   dataI = data_eytar->Integral();
   simcI = simc_eytar->Integral();
   scale = 1./data_eytar->Integral();
   simc_eytar->Scale(scale);
   data_eytar->Scale(scale);
   simc_eytar->Draw("samesE0");
   data_eytar->Draw("sameshistE0");
   leg5->AddEntry(data_eytar,"Data","f");
   leg5->AddEntry(simc_eytar,"SIMC");
   leg5->AddEntry((TObject*)0, Form("Data Integral: %d", dataI), "");
   leg5->AddEntry((TObject*)0, Form("SIMC Integral: %d", simcI), "");
   leg5->Draw();
   
   c1->cd(2);
   dataI = data_exptar->Integral();
   simcI = simc_exptar->Integral();
   scale = 1./dataI;
   simc_exptar->Scale(scale);
   data_exptar->Scale(scale);
   simc_exptar->Draw("samesE0");
   data_exptar->Draw("sameshistE0");
   leg6->AddEntry(data_exptar,"Data", "f");
   leg6->AddEntry(simc_exptar,"SIMC");
   leg6->AddEntry((TObject*)0, Form("Data Integral: %d", dataI), "");
   leg6->AddEntry((TObject*)0, Form("SIMC Integral: %d", simcI), "");
   leg6->Draw();
   
   c1->cd(3);
   dataI = data_eyptar->Integral();
   simcI = simc_eyptar->Integral();
   scale = 1./dataI;
   simc_eyptar->Scale(scale);
   data_eyptar->Scale(scale);
   simc_eyptar->Draw("samesE0");
   data_eyptar->Draw("sameshistE0");
   leg7->AddEntry(data_eyptar,"Data", "f");
   leg7->AddEntry(simc_eyptar,"SIMC");
   leg7->AddEntry((TObject*)0, Form("Data Integral: %d", dataI), "");
   leg7->AddEntry((TObject*)0, Form("SIMC Integral: %d", simcI), "");
   leg7->Draw();
     
   c1->cd(4);
   dataI = data_edelta->Integral();
   simcI = simc_edelta->Integral();
   scale = 1./dataI;
   simc_edelta->Scale(scale);
   data_edelta->Scale(scale);
   simc_edelta->Draw("samesE0");
   data_edelta->Draw("sameshistE0");
   leg8->AddEntry(data_edelta,"Data", "f");
   leg8->AddEntry(simc_edelta,"SIMC");
   leg8->AddEntry((TObject*)0, Form("Data Integral: %d", dataI), "");
   leg8->AddEntry((TObject*)0, Form("SIMC Integral: %d", simcI), "");
   leg8->Draw();
   
   plots_path = plots_dir + Form("eArm_TargRecon_%d.pdf", run);
   c1->SaveAs(plots_path.c_str());

   if(ratio==1){
   //------------DATA/SIMC RATIO------------
   
   TCanvas *c1_r = new TCanvas("c1_r", "Electron Arm: Target Reconstruction Ratio", 2000, 1000);
   c1_r->Divide(2,2);

   c1_r->cd(1);
   data_eytar->Divide(simc_eytar);
   data_eytar->SetMarkerStyle(20);
   data_eytar->SetMarkerColor(kBlack);
   data_eytar->SetLineWidth(2);
   data_eytar->SetLineColor(kBlack);
   data_eytar->GetYaxis()->SetRangeUser(0.5,1.5);
   data_eytar->Draw();

   c1_r->cd(2);
   data_exptar->Divide(simc_exptar);
   data_exptar->SetMarkerStyle(20);
   data_exptar->SetMarkerColor(kBlack);
   data_exptar->SetLineWidth(2);
   data_exptar->SetLineColor(kBlack);
   data_exptar->GetYaxis()->SetRangeUser(0.5,1.5);
   data_exptar->Draw();
   
   c1_r->cd(3);
   data_eyptar->Divide(simc_eyptar);
   data_eyptar->SetMarkerStyle(20);
   data_eyptar->SetMarkerColor(kBlack);
   data_eyptar->SetLineWidth(2);
   data_eyptar->SetLineColor(kBlack);
   data_eyptar->GetYaxis()->SetRangeUser(0.5,1.5);
   data_eyptar->Draw();

   c1_r->cd(4);
   data_edelta->Divide(simc_edelta);
   data_edelta->SetMarkerStyle(20);
   data_edelta->SetMarkerColor(kBlack);
   data_edelta->SetLineWidth(2);
   data_edelta->SetLineColor(kBlack);
   data_edelta->GetYaxis()->SetRangeUser(0.5,1.5);
   data_edelta->Draw();
   
   plots_path = plots_dir + Form("RATIO_eArm_TargRecon_%d.pdf", run);
   c1_r->SaveAs(plots_path.c_str());
      
   }
   
   //------------------------------------------------------------------------------
   
   
   //-----------------PLOT FOCAL PLANE  Variables SIMC/Data comparison-----------------------

  //Set Legend
   auto leg9 = new TLegend(0.1,0.8,0.28,0.9);
   auto leg10 = new TLegend(0.1,0.8,0.28,0.9);
   auto leg11 = new TLegend(0.1,0.8,0.28,0.9);
   auto leg12 = new TLegend(0.1,0.8,0.28,0.9);
   
   TCanvas *c2 = new TCanvas("c2", "Electron Arm: Focal Plane", 2000, 1000);
   c2->Divide(2,2);

   c2->cd(1);
   dataI = data_exfp->Integral();
   simcI = simc_exfp->Integral();
   scale = 1./data_exfp->Integral();
   simc_exfp->Scale(scale);
   data_exfp->Scale(scale);
   simc_exfp->Draw("samesE0");
   data_exfp->Draw("sameshistE0");
   leg9->AddEntry(data_exfp,"Data","f");
   leg9->AddEntry(simc_exfp,"SIMC");
   leg9->AddEntry((TObject*)0, Form("Data Integral: %d", dataI), "");
   leg9->AddEntry((TObject*)0, Form("SIMC Integral: %d", simcI), "");
   leg9->Draw();
   
   c2->cd(2);
   dataI = data_eyfp->Integral();
   simcI = simc_eyfp->Integral();
   scale = 1./data_eyfp->Integral();
   simc_eyfp->Scale(scale);
   data_eyfp->Scale(scale);
   simc_eyfp->Draw("samesE0");
   data_eyfp->Draw("sameshistE0");
   leg10->AddEntry(data_eyfp,"Data", "f");
   leg10->AddEntry(simc_eyfp,"SIMC");
   leg10->AddEntry((TObject*)0, Form("Data Integral: %d", dataI), "");
   leg10->AddEntry((TObject*)0, Form("SIMC Integral: %d", simcI), "");
   leg10->Draw();

   c2->cd(3);
   dataI = data_expfp->Integral();
   simcI = simc_expfp->Integral();
   scale = 1./data_expfp->Integral();
   simc_expfp->Scale(scale);
   data_expfp->Scale(scale);
   simc_expfp->Draw("samesE0");
   data_expfp->Draw("sameshistE0");
   leg11->AddEntry(data_expfp,"Data", "f");
   leg11->AddEntry(simc_expfp,"SIMC");
   leg11->AddEntry((TObject*)0, Form("Data Integral: %d", dataI), "");
   leg11->AddEntry((TObject*)0, Form("SIMC Integral: %d", simcI), "");
   leg11->Draw();
     
   c2->cd(4);
   dataI = data_eypfp->Integral();
   simcI = simc_eypfp->Integral();
   scale = 1./data_eypfp->Integral();
   simc_eypfp->Scale(scale);
   data_eypfp->Scale(scale);
   simc_eypfp->Draw("samesE0");
   data_eypfp->Draw("sameshistE0");
   leg12->AddEntry(data_eypfp,"Data", "f");
   leg12->AddEntry(simc_eypfp,"SIMC");
   leg12->AddEntry((TObject*)0, Form("Data Integral: %d", dataI), "");
   leg12->AddEntry((TObject*)0, Form("SIMC Integral: %d", simcI), "");
   leg12->Draw();

   plots_path = plots_dir + Form("eArm_FocalPlane_%d.pdf", run);
   c2->SaveAs(plots_path.c_str());                                                                                   

   if(ratio==1){
   TCanvas *c2_r = new TCanvas("c2_r", "Electron Arm: Target Reconstruction Ratio", 2000, 1000);
   c2_r->Divide(2,2);

   c2_r->cd(1);
   data_exfp->Divide(simc_exfp);
   data_exfp->SetMarkerStyle(20);
   data_exfp->SetMarkerColor(kBlack);
   data_exfp->SetLineWidth(2);
   data_exfp->SetLineColor(kBlack);
   data_exfp->GetYaxis()->SetRangeUser(0.5,1.5);
   data_exfp->Draw();

   c2_r->cd(2);
   data_eyfp->Divide(simc_eyfp);
   data_eyfp->SetMarkerStyle(20);
   data_eyfp->SetMarkerColor(kBlack);
   data_eyfp->SetLineWidth(2);
   data_eyfp->SetLineColor(kBlack);
   data_eyfp->GetYaxis()->SetRangeUser(0.5,1.5);
   data_eyfp->Draw();

   c2_r->cd(3);
   data_expfp->Divide(simc_exfp);
   data_expfp->SetMarkerStyle(20);
   data_expfp->SetMarkerColor(kBlack);
   data_expfp->SetLineWidth(2);
   data_expfp->SetLineColor(kBlack);
   data_expfp->GetYaxis()->SetRangeUser(0.5,1.5);
   data_expfp->Draw();

   c2_r->cd(4);
   data_eypfp->Divide(simc_exfp);
   data_eypfp->SetMarkerStyle(20);
   data_eypfp->SetMarkerColor(kBlack);
   data_eypfp->SetLineWidth(2);
   data_eypfp->SetLineColor(kBlack);
   data_eypfp->GetYaxis()->SetRangeUser(0.5,1.5);
   data_eypfp->Draw();

   }
   //----------------------------------------------------------- 
 
   
   //-----------------PLOT KINEMATICS SIMC/Data comparison---------------

//Set Legend

   auto leg13 = new TLegend(0.1,0.8,0.28,0.9);
   auto leg14 = new TLegend(0.1,0.8,0.28,0.9);
   auto leg15 = new TLegend(0.1,0.8,0.28,0.9);
   auto leg16 = new TLegend(0.1,0.8,0.28,0.9);
   auto leg17 = new TLegend(0.1,0.8,0.28,0.9);
   auto leg18 = new TLegend(0.1,0.8,0.28,0.9);
   auto leg19 = new TLegend(0.1,0.8,0.28,0.9);
   auto leg20 = new TLegend(0.1,0.8,0.28,0.9);

   
   TCanvas *c3 = new TCanvas("c3", "Kinematics", 2000, 1000);
   c3->Divide(4,2);
   
   c3->cd(1);
   dataI = data_Q2->Integral();
   simcI = simc_Q2->Integral();
   scale = 1./data_Q2->Integral();
   simc_Q2->Scale(scale);
   data_Q2->Scale(scale);
   simc_Q2->Draw("samesE0");
   data_Q2->Draw("sameshistE0");
   leg13->AddEntry(data_Q2,"Data", "f");
   leg13->AddEntry(simc_Q2,"SIMC");
   leg13->AddEntry((TObject*)0, Form("Data Integral: %d", dataI), "");
   leg13->AddEntry((TObject*)0, Form("SIMC Integral: %d", simcI), "");
   leg13->Draw();
     
   c3->cd(2);
   dataI = data_omega->Integral();
   simcI = simc_omega->Integral();
   scale = 1./data_omega->Integral();
   simc_omega->Scale(scale);
   data_omega->Scale(scale);
   simc_omega->Draw("samesE0");
   data_omega->Draw("sameshistE0");
   leg14->AddEntry(data_omega,"Data", "f");
   leg14->AddEntry(simc_omega,"SIMC");
   leg14->AddEntry((TObject*)0, Form("Data Integral: %d", dataI), "");
   leg14->AddEntry((TObject*)0, Form("SIMC Integral: %d", simcI), "");
   leg14->Draw();

   c3->cd(3);
   dataI = data_W->Integral();
   simcI = simc_W->Integral();
   scale = 1./data_W->Integral();
   simc_W->Scale(scale);
   data_W->Scale(scale);
   data_W->GetXaxis()->SetTitle("Invariant Mass, W [GeV]");
   data_W->GetXaxis()->CenterTitle();
   simc_W->Draw("samesE0");
   data_W->Draw("sameshistE0");
   leg15->AddEntry(data_W,"Data", "f");
   leg15->AddEntry(simc_W,"SIMC");
   leg15->AddEntry((TObject*)0, Form("Data Integral: %d", dataI), "");
   leg15->AddEntry((TObject*)0, Form("SIMC Integral: %d", simcI), "");
   leg15->Draw();

   c3->cd(4);
   dataI = data_thq->Integral();
   simcI = simc_thq->Integral();
   scale = 1./data_thq->Integral();
   simc_thq->Scale(scale);
   data_thq->Scale(scale);
   data_thq->GetXaxis()->SetTitle("q-vector Angle, #theta_{q} [deg]");
   data_thq->GetXaxis()->CenterTitle();
   simc_thq->Draw("samesE0");
   data_thq->Draw("sameshistE0");
   leg16->AddEntry(data_thq,"Data", "f");
   leg16->AddEntry(simc_thq,"SIMC");
   leg16->AddEntry((TObject*)0, Form("Data Integral: %d", dataI), "");
   leg16->AddEntry((TObject*)0, Form("SIMC Integral: %d", simcI), "");
   leg16->Draw();

   c3->cd(5);
   dataI = data_xbj->Integral();
   simcI = simc_xbj->Integral();
   scale = 1./data_xbj->Integral();
   simc_xbj->Scale(scale);
   data_xbj->Scale(scale);
   simc_xbj->Draw("samesE0");
   data_xbj->Draw("sameshistE0");
   leg17->AddEntry(data_xbj,"Data","f");
   leg17->AddEntry(simc_xbj,"SIMC");
   leg17->AddEntry((TObject*)0, Form("Data Integral: %d", dataI), "");
   leg17->AddEntry((TObject*)0, Form("SIMC Integral: %d", simcI), "");
   leg17->Draw();

   c3->cd(6);
   dataI = data_th_elec->Integral();
   simcI = simc_th_elec->Integral();
   scale = 1./data_th_elec->Integral();
   simc_th_elec->Scale(scale);
   data_th_elec->Scale(scale);
   data_th_elec->GetXaxis()->SetTitle("Electron Scatt. Angle, #theta_{e} [deg]");
   data_th_elec->GetXaxis()->CenterTitle();
   simc_th_elec->Draw("samesE0");
   data_th_elec->Draw("sameshistE0");
   leg18->AddEntry(data_th_elec,"Data","f");
   leg18->AddEntry(simc_th_elec,"SIMC");
   leg18->AddEntry((TObject*)0, Form("Data Integral: %d", dataI), "");
   leg18->AddEntry((TObject*)0, Form("SIMC Integral: %d", simcI), "");
   leg18->Draw();

   c3->cd(7);
   dataI = data_kf->Integral();
   simcI = simc_kf->Integral();
   scale = 1./data_kf->Integral();
   simc_kf->Scale(scale);
   data_kf->Scale(scale);
   data_kf->GetXaxis()->SetTitle("Electron Final Momentum, k_{f} [GeV/c] ");
   data_kf->GetXaxis()->CenterTitle();   
   simc_kf->Draw("samesE0");
   data_kf->Draw("sameshistE0");
   leg19->AddEntry(data_kf,"Data","f");
   leg19->AddEntry(simc_kf,"SIMC");
   leg19->AddEntry((TObject*)0, Form("Data Integral: %d", dataI), "");
   leg19->AddEntry((TObject*)0, Form("SIMC Integral: %d", simcI), "");
   leg19->Draw();

   c3->cd(8);
   dataI = data_emiss->Integral();
   simcI = simc_emiss->Integral();
   scale = 1./data_emiss->Integral();
   simc_emiss->Scale(scale);
   data_emiss->Scale(scale);
   data_emiss->GetXaxis()->SetTitle("Missing Energy, E_{m} [GeV/c] ");
   data_emiss->GetXaxis()->CenterTitle();   
   simc_emiss->Draw("samesE0");
   data_emiss->Draw("sameshistE0");
   leg20->AddEntry(data_emiss,"Data","f");
   leg20->AddEntry(simc_emiss,"SIMC");
   leg20->AddEntry((TObject*)0, Form("Data Integral: %d", dataI), "");
   leg20->AddEntry((TObject*)0, Form("SIMC Integral: %d", simcI), "");
   leg20->Draw();

   plots_path = plots_dir + Form("Kinematics1_%d.pdf", run);
   c3->SaveAs(plots_path.c_str());                                                                   


   if(ratio==1)
     {
   //----------DATA/SIMC RATIO------
   TCanvas *c3_r = new TCanvas("c3_r", "Kinematics", 2000, 1000);
   c3_r->Divide(4,2);

   c3_r->cd(1);
   data_Q2->Divide(simc_Q2);
   data_Q2->SetMarkerStyle(20);
   data_Q2->SetMarkerColor(kBlack);
   data_Q2->SetLineWidth(2);
   data_Q2->SetLineColor(kBlack);
   data_Q2->GetYaxis()->SetRangeUser(0.5,1.5);
   data_Q2->Draw();

   c3_r->cd(2);
   data_omega->Divide(simc_omega);
   data_omega->SetMarkerStyle(20);
   data_omega->SetMarkerColor(kBlack);
   data_omega->SetLineWidth(2);
   data_omega->SetLineColor(kBlack);
   data_omega->GetYaxis()->SetRangeUser(0.5,1.5);
   data_omega->Draw();

   c3_r->cd(3);
   data_W->Divide(simc_W);
   data_W->SetMarkerStyle(20);
   data_W->SetMarkerColor(kBlack);
   data_W->SetLineWidth(2);
   data_W->SetLineColor(kBlack);
   data_W->GetYaxis()->SetRangeUser(0.5,1.5);
   data_W->Draw();

   c3_r->cd(4);
   data_thq->Divide(simc_thq);
   data_thq->SetMarkerStyle(20);
   data_thq->SetMarkerColor(kBlack);
   data_thq->SetLineWidth(2);
   data_thq->SetLineColor(kBlack);
   data_thq->GetYaxis()->SetRangeUser(0.5,1.5);
   data_thq->Draw();

   c3_r->cd(5);
   data_xbj->Divide(simc_xbj);
   data_xbj->SetMarkerStyle(20);
   data_xbj->SetMarkerColor(kBlack);
   data_xbj->SetLineWidth(2);
   data_xbj->SetLineColor(kBlack);
   data_xbj->GetYaxis()->SetRangeUser(0.5,1.5);
   data_xbj->Draw();

   c3_r->cd(6);
   data_th_elec->Divide(simc_th_elec);
   data_th_elec->SetMarkerStyle(20);
   data_th_elec->SetMarkerColor(kBlack);
   data_th_elec->SetLineWidth(2);
   data_th_elec->SetLineColor(kBlack);
   data_th_elec->GetYaxis()->SetRangeUser(0.5,1.5);
   data_th_elec->Draw();

   c3_r->cd(7);
   data_kf->Divide(simc_kf);
   data_kf->SetMarkerStyle(20);
   data_kf->SetMarkerColor(kBlack);
   data_kf->SetLineWidth(2);
   data_kf->SetLineColor(kBlack);
   data_kf->GetYaxis()->SetRangeUser(0.5,1.5);
   data_kf->Draw();

   c3_r->cd(8);
   data_emiss->Divide(simc_emiss);
   data_emiss->SetMarkerStyle(20);
   data_emiss->SetMarkerColor(kBlack);
   data_emiss->SetLineWidth(2);
   data_emiss->SetLineColor(kBlack);
   data_emiss->GetYaxis()->SetRangeUser(0.5,1.5);
   data_emiss->Draw();

     }
   
   //---------------------------
   //Plot Additional Kinematics
   
   auto leg_Pm = new TLegend(0.1,0.8,0.28,0.9);
   auto leg_Pf = new TLegend(0.1,0.8,0.28,0.9);
   auto leg_thp = new TLegend(0.1,0.8,0.28,0.9);
   auto leg_q = new TLegend(0.1,0.8,0.28,0.9);
   auto leg_thpq = new TLegend(0.1,0.8,0.28,0.9);
   auto leg_Pmx = new TLegend(0.1,0.8,0.28,0.9);
   auto leg_Pmy = new TLegend(0.1,0.8,0.28,0.9);
   auto leg_Pmz = new TLegend(0.1,0.8,0.28,0.9);


   //Create A Canvas to store kinematic variable comparisons
   TCanvas *ck2 = new TCanvas("ck2", "Kinematics-2", 2000, 1000);
   ck2->Divide(4,2);
   
   ck2->cd(1);
   dataI = data_Pm->Integral();
   simcI = simc_Pm->Integral();
   scale = 1./data_Pm->Integral();
   simc_Pm->Scale(scale);
   data_Pm->Scale(scale);
   data_Pm->GetXaxis()->SetTitle("Missing Momentum, P_{miss} [GeV]");
   data_Pm->GetXaxis()->CenterTitle();
   simc_Pm->Draw("samesE0");
   data_Pm->Draw("sameshistE0");
   leg_Pm->AddEntry(data_Pm,"Data", "f");
   leg_Pm->AddEntry(simc_Pm,"SIMC");
   leg_Pm->AddEntry((TObject*)0, Form("Data Integral: %d", dataI), "");
   leg_Pm->AddEntry((TObject*)0, Form("SIMC Integral: %d", simcI), "");
   leg_Pm->Draw();

   ck2->cd(2);
   dataI = data_Pf->Integral();
   simcI = simc_Pf->Integral();
   scale = 1./data_Pf->Integral();
   simc_Pf->Scale(scale);
   data_Pf->Scale(scale);
   data_Pf->GetXaxis()->SetTitle("Proton Momentum, P_{p} [GeV]");
   data_Pf->GetXaxis()->CenterTitle();
   simc_Pf->Draw("samesE0");
   data_Pf->Draw("sameshistE0");
   leg_Pf->AddEntry(data_Pf,"Data", "f");
   leg_Pf->AddEntry(simc_Pf,"SIMC");
   leg_Pf->AddEntry((TObject*)0, Form("Data Integral: %d", dataI), "");
   leg_Pf->AddEntry((TObject*)0, Form("SIMC Integral: %d", simcI), "");
   leg_Pf->Draw();

   ck2->cd(3);
   dataI = data_th_prot->Integral();
   simcI = simc_th_prot->Integral();
   scale = 1./data_th_prot->Integral();
   simc_th_prot->Scale(scale);
   data_th_prot->Scale(scale);
   data_th_prot->GetXaxis()->SetTitle("Proton Scatt. Angle, #theta_{p} [deg]");
   data_th_prot->GetXaxis()->CenterTitle();
   simc_th_prot->Draw("samesE0");
   data_th_prot->Draw("sameshistE0");
   leg_thp->AddEntry(data_th_prot,"Data", "f");
   leg_thp->AddEntry(simc_th_prot,"SIMC");
   leg_thp->AddEntry((TObject*)0, Form("Data Integral: %d", dataI), "");
   leg_thp->AddEntry((TObject*)0, Form("SIMC Integral: %d", simcI), "");
   leg_thp->Draw();

   ck2->cd(4);
   dataI = data_q->Integral();
   simcI = simc_q->Integral();
   scale = 1./data_q->Integral();
   simc_q->Scale(scale);
   data_q->Scale(scale);
   data_q->GetXaxis()->SetTitle("q-Vector Magnitude, |q| [GeV]");
   data_q->GetXaxis()->CenterTitle();
   simc_q->Draw("samesE0");
   data_q->Draw("sameshistE0");
   leg_q->AddEntry(data_q,"Data", "f");
   leg_q->AddEntry(simc_q,"SIMC");
   leg_q->AddEntry((TObject*)0, Form("Data Integral: %d", dataI), "");
   leg_q->AddEntry((TObject*)0, Form("SIMC Integral: %d", simcI), "");
   leg_q->Draw();


   ck2->cd(5);
   dataI = data_thpq->Integral();
   simcI = simc_thpq->Integral();
   scale = 1./data_thpq->Integral();
   simc_thpq->Scale(scale);
   data_thpq->Scale(scale);
   data_thpq->GetXaxis()->SetTitle("Proton-qVec. Angle, #theta_{pq} [deg]");
   data_thpq->GetXaxis()->CenterTitle();
   simc_thpq->Draw("samesE0");
   data_thpq->Draw("sameshistE0");
   leg_thpq->AddEntry(data_thpq,"Data", "f");
   leg_thpq->AddEntry(simc_thpq,"SIMC");
   leg_thpq->AddEntry((TObject*)0, Form("Data Integral: %d", dataI), "");
   leg_thpq->AddEntry((TObject*)0, Form("SIMC Integral: %d", simcI), "");
   leg_thpq->Draw();
  
   ck2->cd(6);
   dataI = data_Pmx->Integral();
   simcI = simc_Pmx->Integral();
   scale = 1./data_Pmx->Integral();
   simc_Pmx->Scale(scale);
   data_Pmx->Scale(scale);
   data_Pmx->GetXaxis()->SetTitle("Missing Momentum X-comp., Pm_{x} [GeV]");
   data_Pmx->GetXaxis()->CenterTitle();
   simc_Pmx->Draw("samesE0");
   data_Pmx->Draw("sameshistE0");
   leg_Pmx->AddEntry(data_Pmx,"Data", "f");
   leg_Pmx->AddEntry(simc_Pmx,"SIMC");
   leg_Pmx->AddEntry((TObject*)0, Form("Data Integral: %d", dataI), "");
   leg_Pmx->AddEntry((TObject*)0, Form("SIMC Integral: %d", simcI), "");
   leg_Pmx->Draw();

   ck2->cd(7);
   dataI = data_Pmy->Integral();
   simcI = simc_Pmy->Integral();
   scale = 1./data_Pmy->Integral();
   simc_Pmy->Scale(scale);
   data_Pmy->Scale(scale);
   data_Pmy->GetXaxis()->SetTitle("Missing Momentum Y-comp., Pm_{y} [GeV]");
   data_Pmy->GetXaxis()->CenterTitle();
   simc_Pmy->Draw("samesE0");
   data_Pmy->Draw("sameshistE0");
   leg_Pmy->AddEntry(data_Pmy,"Data", "f");
   leg_Pmy->AddEntry(simc_Pmy,"SIMC");
   leg_Pmy->AddEntry((TObject*)0, Form("Data Integral: %d", dataI), "");
   leg_Pmy->AddEntry((TObject*)0, Form("SIMC Integral: %d", simcI), "");
   leg_Pmy->Draw();

   ck2->cd(8);
   dataI = data_Pmz->Integral();
   simcI = simc_Pmz->Integral();
   scale = 1./data_Pmz->Integral();
   simc_Pmz->Scale(scale);
   data_Pmz->Scale(scale);
   data_Pmz->GetXaxis()->SetTitle("Missing Momentum Z-comp., Pm_{z} [GeV]");
   data_Pmz->GetXaxis()->CenterTitle();
   simc_Pmz->Draw("samesE0");
   data_Pmz->Draw("sameshistE0");
   leg_Pmz->AddEntry(data_Pmz,"Data", "f");
   leg_Pmz->AddEntry(simc_Pmz,"SIMC");
   leg_Pmz->AddEntry((TObject*)0, Form("Data Integral: %d", dataI), "");
   leg_Pmz->AddEntry((TObject*)0, Form("SIMC Integral: %d", simcI), "");
   leg_Pmz->Draw();
   
   plots_path = plots_dir + Form("Kinematics2_%d.pdf", run);
   ck2->SaveAs(plots_path.c_str());                                                                   

   if(ratio==1){
   //----------DATA/SIMC RATIO------
   TCanvas *ck2_r = new TCanvas("ck2_r", "Kinematics-2 Ratio", 2000, 1000);
   ck2_r->Divide(4,2);

   ck2_r->cd(1);
   data_Pm->Divide(simc_Pm);
   data_Pm->SetMarkerStyle(20);
   data_Pm->SetMarkerColor(kBlack);
   data_Pm->SetLineWidth(2);
   data_Pm->SetLineColor(kBlack);
   data_Pm->GetYaxis()->SetRangeUser(0.5,1.5);
   data_Pm->Draw();

   ck2_r->cd(2);
   data_Pf->Divide(simc_Pf);
   data_Pf->SetMarkerStyle(20);
   data_Pf->SetMarkerColor(kBlack);
   data_Pf->SetLineWidth(2);
   data_Pf->SetLineColor(kBlack);
   data_Pf->GetYaxis()->SetRangeUser(0.5,1.5);
   data_Pf->Draw();
   
   ck2_r->cd(3);
   data_th_prot->Divide(simc_th_prot);
   data_th_prot->SetMarkerStyle(20);
   data_th_prot->SetMarkerColor(kBlack);
   data_th_prot->SetLineWidth(2);
   data_th_prot->SetLineColor(kBlack);
   data_th_prot->GetYaxis()->SetRangeUser(0.5,1.5);
   data_th_prot->Draw();

   ck2_r->cd(4);
   data_q->Divide(simc_q);
   data_q->SetMarkerStyle(20);
   data_q->SetMarkerColor(kBlack);
   data_q->SetLineWidth(2);
   data_q->SetLineColor(kBlack);
   data_q->GetYaxis()->SetRangeUser(0.5,1.5);
   data_q->Draw();

   ck2_r->cd(5);
   data_thpq->Divide(simc_thpq);
   data_thpq->SetMarkerStyle(20);
   data_thpq->SetMarkerColor(kBlack);
   data_thpq->SetLineWidth(2);
   data_thpq->SetLineColor(kBlack);
   data_thpq->GetYaxis()->SetRangeUser(0.5,1.5);
   data_thpq->Draw();

   ck2_r->cd(6);
   data_Pmx->Divide(simc_Pmx);
   data_Pmx->SetMarkerStyle(20);
   data_Pmx->SetMarkerColor(kBlack);
   data_Pmx->SetLineWidth(2);
   data_Pmx->SetLineColor(kBlack);
   data_Pmx->GetYaxis()->SetRangeUser(0.5,1.5);
   data_Pmx->Draw();

   ck2_r->cd(7);
   data_Pmy->Divide(simc_Pmy);
   data_Pmy->SetMarkerStyle(20);
   data_Pmy->SetMarkerColor(kBlack);
   data_Pmy->SetLineWidth(2);
   data_Pmy->SetLineColor(kBlack);
   data_Pmy->GetYaxis()->SetRangeUser(0.5,1.5);
   data_Pmy->Draw();

   ck2_r->cd(8);
   data_Pmz->Divide(simc_Pmz);
   data_Pmz->SetMarkerStyle(20);
   data_Pmz->SetMarkerColor(kBlack);
   data_Pmz->SetLineWidth(2);
   data_Pmz->SetLineColor(kBlack);
   data_Pmz->GetYaxis()->SetRangeUser(0.5,1.5);
   data_Pmz->Draw();

   }
   
   //-----------------PLOT TARGET  Variables SIMC/Data comparison-----------------------

  //Set Legend
   auto leghxt = new TLegend(0.1,0.8,0.28,0.9);                          
   auto leghyt = new TLegend(0.1,0.8,0.28,0.9);  
   auto leghzt = new TLegend(0.1,0.8,0.28,0.9);                                                          
                                                                                                                                                          
   auto legpxt = new TLegend(0.1,0.8,0.28,0.9);                              
   auto legpyt = new TLegend(0.1,0.8,0.28,0.9);                                                                       
   auto legpzt = new TLegend(0.1,0.8,0.28,0.9);   

   TCanvas *c4a = new TCanvas("c4a", "HMS Target Variables", 2000, 1000);
   c4a->Divide(3,1);

   c4a->cd(1);
   dataI = data_xtarH->Integral();
   simcI = simc_xtarH->Integral();
   scale = 1./data_xtarH->Integral();
   simc_xtarH->Scale(scale);
   data_xtarH->Scale(scale);
   simc_xtarH->Draw("samesE0");
   data_xtarH->Draw("sameshistE0");
   leghxt->AddEntry(data_xtarH,"Data","f");
   leghxt->AddEntry(simc_xtarH,"SIMC");
   leghxt->AddEntry((TObject*)0, Form("Data Integral: %d", dataI), "");
   leghxt->AddEntry((TObject*)0, Form("SIMC Integral: %d", simcI), "");
   leghxt->Draw();
  
   c4a->cd(2);
   dataI = data_ytarH->Integral();
   simcI = simc_ytarH->Integral();
   scale = 1./data_ytarH->Integral();
   simc_ytarH->Scale(scale);
   data_ytarH->Scale(scale);
   simc_ytarH->Draw("samesE0");
   data_ytarH->Draw("sameshistE0");
   leghyt->AddEntry(data_ytarH,"Data","f");
   leghyt->AddEntry(simc_ytarH,"SIMC");
   leghyt->AddEntry((TObject*)0, Form("Data Integral: %d", dataI), "");
   leghyt->AddEntry((TObject*)0, Form("SIMC Integral: %d", simcI), "");
   leghyt->Draw();

   c4a->cd(3);
   dataI = data_ztarH->Integral();
   simcI = simc_ztarH->Integral();
   scale = 1./data_ztarH->Integral();
   simc_ztarH->Scale(scale);
   data_ztarH->Scale(scale);
   simc_ztarH->Draw("samesE0");
   data_ztarH->Draw("sameshistE0");
   leghzt->AddEntry(data_ztarH,"Data","f");
   leghzt->AddEntry(simc_ztarH,"SIMC");
   leghzt->AddEntry((TObject*)0, Form("Data Integral: %d", dataI), "");
   leghzt->AddEntry((TObject*)0, Form("SIMC Integral: %d", simcI), "");
   leghzt->Draw();
   
   plots_path = plots_dir + Form("hArm_TargVar_%d.pdf", run);
   c4a->SaveAs(plots_path.c_str());                                                                                              

   if(ratio==1){
   //-------DATA/SIMC RATIO-----
   TCanvas *c4a_r = new TCanvas("c4a_r", "HMS Target Variables", 2000, 2000);
   c4a_r->Divide(3,1);

   c4a_r->cd(1);
   data_xtarH->Divide(simc_xtarH);
   data_xtarH->SetMarkerStyle(20);
   data_xtarH->SetMarkerColor(kBlack);
   data_xtarH->SetLineWidth(2);
   data_xtarH->SetLineColor(kBlack);
   data_xtarH->GetYaxis()->SetRangeUser(0.5,1.5);
   data_xtarH->Draw();

   c4a_r->cd(2);
   data_ytarH->Divide(simc_ytarH);
   data_ytarH->SetMarkerStyle(20);
   data_ytarH->SetMarkerColor(kBlack);
   data_ytarH->SetLineWidth(2);
   data_ytarH->SetLineColor(kBlack);
   data_ytarH->GetYaxis()->SetRangeUser(0.5,1.5);
   data_ytarH->Draw();

   c4a_r->cd(3);
   data_ztarH->Divide(simc_ztarH);
   data_ztarH->SetMarkerStyle(20);
   data_ztarH->SetMarkerColor(kBlack);
   data_ztarH->SetLineWidth(2);
   data_ztarH->SetLineColor(kBlack);
   data_ztarH->GetYaxis()->SetRangeUser(0.5,1.5);
   data_ztarH->Draw();

   }
   
   //------------SHMS TARGET VARIABLES-------
   
   TCanvas *c4b = new TCanvas("c4b", "SHMS Target Variables", 2000, 1000);
   c4b->Divide(3,1);

   c4b->cd(1);
   dataI = data_xtarP->Integral();
   simcI = simc_xtarP->Integral();
   scale = 1./data_xtarP->Integral();
   simc_xtarP->Scale(scale);
   data_xtarP->Scale(scale);
   simc_xtarP->Draw("samesE0");
   data_xtarP->Draw("sameshistE0E0");
   legpxt->AddEntry(data_xtarP,"Data","f");
   legpxt->AddEntry(simc_xtarP,"SIMC");
   legpxt->AddEntry((TObject*)0, Form("Data Integral: %d", dataI), "");
   legpxt->AddEntry((TObject*)0, Form("SIMC Integral: %d", simcI), "");
   legpxt->Draw();
  
   c4b->cd(2);
   dataI = data_ytarP->Integral();
   simcI = simc_ytarP->Integral();
   scale = 1./data_ytarP->Integral();
   simc_ytarP->Scale(scale);
   data_ytarP->Scale(scale);
   simc_ytarP->Draw("samesE0");
   data_ytarP->Draw("sameshistE0E0");
   legpyt->AddEntry(data_ytarP,"Data","f");
   legpyt->AddEntry(simc_ytarP,"SIMC");
   legpyt->AddEntry((TObject*)0, Form("Data Integral: %d", dataI), "");
   legpyt->AddEntry((TObject*)0, Form("SIMC Integral: %d", simcI), "");
   legpyt->Draw();

   c4b->cd(3);
   dataI = data_ztarP->Integral();
   simcI = simc_ztarP->Integral();
   scale = 1./data_ztarP->Integral();
   simc_ztarP->Scale(scale);
   data_ztarP->Scale(scale);
   simc_ztarP->Draw("samesE0");
   data_ztarP->Draw("sameshistE0E0");
   legpzt->AddEntry(data_ztarP,"Data","f");
   legpzt->AddEntry(simc_ztarP,"SIMC");
   legpzt->AddEntry((TObject*)0, Form("Data Integral: %d", dataI), "");
   legpzt->AddEntry((TObject*)0, Form("SIMC Integral: %d", simcI), "");
   legpzt->Draw();
   
   plots_path = plots_dir + Form("pArm_TargVar_%d.pdf", run);
   c4b->SaveAs(plots_path.c_str());      

   if(ratio==1){
 //-------DATA/SIMC RATIO-----
   TCanvas *c4b_r = new TCanvas("c4b_r", "RATIO SHMS Target Variables", 2000, 2000);
   c4b_r->Divide(3,1);

   c4b_r->cd(1);
   data_xtarP->Divide(simc_xtarP);
   data_xtarP->SetMarkerStyle(20);
   data_xtarP->SetMarkerColor(kBlack);
   data_xtarP->SetLineWidth(2);
   data_xtarP->SetLineColor(kBlack);
   data_xtarP->GetYaxis()->SetRangeUser(0.5,1.5);
   data_xtarP->Draw();

   c4b_r->cd(2);
   data_ytarP->Divide(simc_ytarP);
   data_ytarP->SetMarkerStyle(20);
   data_ytarP->SetMarkerColor(kBlack);
   data_ytarP->SetLineWidth(2);
   data_ytarP->SetLineColor(kBlack);
   data_ytarP->GetYaxis()->SetRangeUser(0.5,1.5);
   data_ytarP->Draw();

   c4b_r->cd(3);
   data_ztarP->Divide(simc_ztarP);
   data_ztarP->SetMarkerStyle(20);
   data_ztarP->SetMarkerColor(kBlack);
   data_ztarP->SetLineWidth(2);
   data_ztarP->SetLineColor(kBlack);
   data_ztarP->GetYaxis()->SetRangeUser(0.5,1.5);
   data_ztarP->Draw();
   }
   
   //-----------------PLOT Target Reconstructed Variables SIMC/Data comparison-----------------------
 
   //Set Legend
   auto htr_l1 = new TLegend(0.1,0.8,0.28,0.9);
   auto htr_l2 = new TLegend(0.1,0.8,0.28,0.9);
   auto htr_l3 = new TLegend(0.1,0.8,0.28,0.9);
   auto htr_l4 = new TLegend(0.1,0.8,0.28,0.9);
   
   //Create A Canvas to store Target Recon. variable comparisons in HADRON ARM
   
   TCanvas *htr = new TCanvas("htr", "Hadron Arm: Target Reconstruction", 2000, 1000);
   htr->Divide(2,2);

   htr->cd(1);
   dataI = data_hytar->Integral();
   simcI = simc_hytar->Integral();
   scale = 1./data_hytar->Integral();
   simc_hytar->Scale(scale);
   data_hytar->Scale(scale);
   simc_hytar->Draw("samesE0");
   data_hytar->Draw("sameshistE0");
   htr_l1->AddEntry(data_hytar,"Data","f");
   htr_l1->AddEntry(simc_hytar,"SIMC");
   htr_l1->AddEntry((TObject*)0, Form("Data Integral: %d", dataI), "");
   htr_l1->AddEntry((TObject*)0, Form("SIMC Integral: %d", simcI), "");
   htr_l1->Draw();

   htr->cd(2);
   dataI = data_hxptar->Integral();
   simcI = simc_hxptar->Integral();
   scale = 1./dataI;
   simc_hxptar->Scale(scale);
   data_hxptar->Scale(scale);
   simc_hxptar->Draw("samesE0");
   data_hxptar->Draw("sameshistE0");
   htr_l2->AddEntry(data_hxptar,"Data", "f");
   htr_l2->AddEntry(simc_hxptar,"SIMC");
   htr_l2->AddEntry((TObject*)0, Form("Data Integral: %d", dataI), "");
   htr_l2->AddEntry((TObject*)0, Form("SIMC Integral: %d", simcI), "");
   htr_l2->Draw();

   htr->cd(3);
   dataI = data_hyptar->Integral();
   simcI = simc_hyptar->Integral();
   scale = 1./dataI;
   simc_hyptar->Scale(scale);
   data_hyptar->Scale(scale);
   simc_hyptar->Draw("samesE0");
   data_hyptar->Draw("sameshistE0");
   htr_l3->AddEntry(data_hyptar,"Data", "f");
   htr_l3->AddEntry(simc_hyptar,"SIMC");
   htr_l3->AddEntry((TObject*)0, Form("Data Integral: %d", dataI), "");
   htr_l3->AddEntry((TObject*)0, Form("SIMC Integral: %d", simcI), "");
   htr_l3->Draw();
     
   htr->cd(4);
   dataI = data_hdelta->Integral();
   simcI = simc_hdelta->Integral();
   scale = 1./dataI;
   simc_hdelta->Scale(scale);
   data_hdelta->Scale(scale);
   simc_hdelta->Draw("samesE0");
   data_hdelta->Draw("sameshistE0");
   htr_l4->AddEntry(data_hdelta,"Data", "f");
   htr_l4->AddEntry(simc_hdelta,"SIMC");
   htr_l4->AddEntry((TObject*)0, Form("Data Integral: %d", dataI), "");
   htr_l4->AddEntry((TObject*)0, Form("SIMC Integral: %d", simcI), "");
   htr_l4->Draw();
   
   
   plots_path = plots_dir + Form("hArm_TargRecon_%d.pdf", run);
   htr->SaveAs(plots_path.c_str());

   if(ratio==1){
    //------------DATA/SIMC RATIO------------
   
   TCanvas *htr_R = new TCanvas("htr_R", "Hadron Arm: Target Reconstruction Ratio", 2000, 1000);
   htr_R->Divide(2,2);

   htr_R->cd(1);
   data_hytar->Divide(simc_hytar);
   data_hytar->SetMarkerStyle(20);
   data_hytar->SetMarkerColor(kBlack);
   data_hytar->SetLineWidth(2);
   data_hytar->SetLineColor(kBlack);
   data_hytar->GetYaxis()->SetRangeUser(0.5,1.5);
   data_hytar->Draw();

   htr_R->cd(2);
   data_hxptar->Divide(simc_hxptar);
   data_hxptar->SetMarkerStyle(20);
   data_hxptar->SetMarkerColor(kBlack);
   data_hxptar->SetLineWidth(2);
   data_hxptar->SetLineColor(kBlack);
   data_hxptar->GetYaxis()->SetRangeUser(0.5,1.5);
   data_hxptar->Draw();
   
   htr_R->cd(3);
   data_hyptar->Divide(simc_hyptar);
   data_hyptar->SetMarkerStyle(20);
   data_hyptar->SetMarkerColor(kBlack);
   data_hyptar->SetLineWidth(2);
   data_hyptar->SetLineColor(kBlack);
   data_hyptar->GetYaxis()->SetRangeUser(0.5,1.5);
   data_hyptar->Draw();

   htr_R->cd(4);
   data_hdelta->Divide(simc_hdelta);
   data_hdelta->SetMarkerStyle(20);
   data_hdelta->SetMarkerColor(kBlack);
   data_hdelta->SetLineWidth(2);
   data_hdelta->SetLineColor(kBlack);
   data_hdelta->GetYaxis()->SetRangeUser(0.5,1.5);
   data_hdelta->Draw();
   }


   
   //------------------------------------------------------------------------------

   
   //-----------------PLOT FOCAL PLANE  Variables SIMC/Data comparison-----------------------

   //Set Legend
   auto hfp_l1 = new TLegend(0.1,0.8,0.28,0.9);
   auto hfp_l2 = new TLegend(0.1,0.8,0.28,0.9);
   auto hfp_l3 = new TLegend(0.1,0.8,0.28,0.9);
   auto hfp_l4 = new TLegend(0.1,0.8,0.28,0.9);

   TCanvas *hfp = new TCanvas("hfp", "Hadron Arm: Focal Plane", 2000, 1000);
   hfp->Divide(2,2);

   hfp->cd(1);
   dataI = data_hxfp->Integral();
   simcI = simc_hxfp->Integral();
   scale = 1./data_hxfp->Integral();
   simc_hxfp->Scale(scale);
   data_hxfp->Scale(scale);
   simc_hxfp->Draw("samesE0");
   data_hxfp->Draw("sameshistE0");
   hfp_l1->AddEntry(data_hxfp,"Data","f");
   hfp_l1->AddEntry(simc_hxfp,"SIMC");
   hfp_l1->AddEntry((TObject*)0, Form("Data Integral: %d", dataI), "");
   hfp_l1->AddEntry((TObject*)0, Form("SIMC Integral: %d", simcI), "");
   hfp_l1->Draw();
   
   hfp->cd(2);
   dataI = data_hyfp->Integral();
   simcI = simc_hyfp->Integral();
   scale = 1./data_hyfp->Integral();
   simc_hyfp->Scale(scale);
   data_hyfp->Scale(scale);
   simc_hyfp->Draw("samesE0");
   data_hyfp->Draw("sameshistE0");
   hfp_l2->AddEntry(data_hyfp,"Data", "f");
   hfp_l2->AddEntry(simc_hyfp,"SIMC");
   hfp_l2->AddEntry((TObject*)0, Form("Data Integral: %d", dataI), "");
   hfp_l2->AddEntry((TObject*)0, Form("SIMC Integral: %d", simcI), "");
   hfp_l2->Draw();

   hfp->cd(3);
   dataI = data_hxpfp->Integral();
   simcI = simc_hxpfp->Integral();
   scale = 1./data_hxpfp->Integral();
   simc_hxpfp->Scale(scale);
   data_hxpfp->Scale(scale);
   simc_hxpfp->Draw("samesE0");
   data_hxpfp->Draw("sameshistE0");
   hfp_l3->AddEntry(data_hxpfp,"Data", "f");
   hfp_l3->AddEntry(simc_hxpfp,"SIMC");
   hfp_l3->AddEntry((TObject*)0, Form("Data Integral: %d", dataI), "");
   hfp_l3->AddEntry((TObject*)0, Form("SIMC Integral: %d", simcI), "");
   hfp_l3->Draw();
     
   hfp->cd(4);
   dataI = data_hypfp->Integral();
   simcI = simc_hypfp->Integral();
   scale = 1./data_hypfp->Integral();
   simc_hypfp->Scale(scale);
   data_hypfp->Scale(scale);
   simc_hypfp->Draw("samesE0");
   data_hypfp->Draw("sameshistE0");
   hfp_l4->AddEntry(data_hypfp,"Data", "f");
   hfp_l4->AddEntry(simc_hypfp,"SIMC");
   hfp_l4->AddEntry((TObject*)0, Form("Data Integral: %d", dataI), "");
   hfp_l4->AddEntry((TObject*)0, Form("SIMC Integral: %d", simcI), "");
   hfp_l4->Draw();
   
   plots_path = plots_dir + Form("hArm_FocalPlane_%d.pdf", run);
   hfp->SaveAs(plots_path.c_str());                                                                                   

   if(ratio==1){
     
   TCanvas *hfp_r = new TCanvas("hfp_r", "Hadron Arm: Focal Plane Ratio", 2000, 1000);
   hfp_r->Divide(2,2);

   hfp_r->cd(1);
   data_hxfp->Divide(simc_hxfp);
   data_hxfp->SetMarkerStyle(20);
   data_hxfp->SetMarkerColor(kBlack);
   data_hxfp->SetLineWidth(2);
   data_hxfp->SetLineColor(kBlack);
   data_hxfp->GetYaxis()->SetRangeUser(0.5,1.5);
   data_hxfp->Draw();

   hfp_r->cd(2);
   data_hyfp->Divide(simc_hyfp);
   data_hyfp->SetMarkerStyle(20);
   data_hyfp->SetMarkerColor(kBlack);
   data_hyfp->SetLineWidth(2);
   data_hyfp->SetLineColor(kBlack);
   data_hyfp->GetYaxis()->SetRangeUser(0.5,1.5);
   data_hyfp->Draw();

   hfp_r->cd(3);
   data_hxpfp->Divide(simc_hxpfp);
   data_hxpfp->SetMarkerStyle(20);
   data_hxpfp->SetMarkerColor(kBlack);
   data_hxpfp->SetLineWidth(2);
   data_hxpfp->SetLineColor(kBlack);
   data_expfp->GetYaxis()->SetRangeUser(0.5,1.5);
   data_hxpfp->Draw();

   hfp_r->cd(4);
   data_hypfp->Divide(simc_hypfp);
   data_hypfp->SetMarkerStyle(20);
   data_hypfp->SetMarkerColor(kBlack);
   data_hypfp->SetLineWidth(2);
   data_hypfp->SetLineColor(kBlack);
   data_hypfp->GetYaxis()->SetRangeUser(0.5,1.5);
   data_hypfp->Draw();


   }
   //----------------------------------------------------------- 
   
   
}
