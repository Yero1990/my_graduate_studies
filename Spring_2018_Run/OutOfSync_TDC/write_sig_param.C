//__________________________________________________________________________
void write_sig_param()
{


  ofstream ofile;
  TString filename = "hdc_sigma_per_wire.param";
  ofile.open(filename);

  TString plane_names[12];
  Int_t nwires[12];

  plane_names[0]="1u1",  nwires[0] = 96;  
  plane_names[1]="1u2",  nwires[1] = 96;
  plane_names[2]="1x1",  nwires[2] = 102;
  plane_names[3]="1x2",  nwires[3] = 102; 
  plane_names[4]="1v2",  nwires[4] = 96;
  plane_names[5]="1v1",  nwires[5] = 96;
  plane_names[6]="2v1",  nwires[6] = 96;
  plane_names[7]="2v2",  nwires[7] = 96;
  plane_names[8]="2x2",  nwires[8] = 102;
  plane_names[9]="2x1",  nwires[9] = 102;
  plane_names[10]="2u2", nwires[10] = 96;
  plane_names[11]="2u1", nwires[11] = 96;
  
  ofile << ";Flag to fix HMS TDC Slot 2 Out-of-Sync problem" << endl;
  ofile << "h_using_sigma_per_wire = 1" << endl;
  ofile << " " << endl;
  
  for (int plane=0; plane<12; plane++) { 
	  
    //write plane headers
    ofile << "hwire_sigma"+plane_names[plane] << "=" << endl;

    Double_t sigma;

    //cout << "Plane: "  << plane << endl;
   for (int wire_num=1; wire_num<=nwires[plane]; wire_num++) 
     {

       if(plane == 0 && wire_num>=81 && wire_num<=96) {sigma = 0.06;}
       else if(plane == 4 && wire_num>=49 && wire_num<=64) {sigma = 0.06;}
       else if(plane == 6 && wire_num>=81 && wire_num<=96) {sigma = 0.06;}
       else if(plane == 10 && wire_num>=49 && wire_num<=64) {sigma = 0.06;}
       else {sigma = 0.02;}
       
       cout << "Plane: " << plane << " ::  wire: " << wire_num << " :: " << "sigma: "<< sigma << endl;
       if (wire_num <= 10) 
	 { 
	   ofile << setprecision(6) << sigma << fixed << ",";
	 }
       else if (wire_num>10 && wire_num <(nwires[plane]))
	 {
	   ofile << setprecision(6) << sigma << ((wire_num+1) % 16 ? ", " : "\n") << fixed;
	 }
       else if (wire_num==nwires[plane]) 
	 {
	   ofile << setprecision(6) << sigma << fixed << endl;
     
	   }
       
     } //END LOOP OVER WIRES

  } //END LOOP OVER PLANES
  
  ofile.close();
  
} //end WriteTZeroParam() Method
