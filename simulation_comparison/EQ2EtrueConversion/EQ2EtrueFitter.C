#include <vector>
#include <iomanip>

void EQ2EtrueFitter(TString geom) {

  gStyle->SetTitleSize(0.08,"t");
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetTitleYSize(0.05);
  //gStyle->SetTitleYOffset(0.1);
  gStyle->SetTitleXSize(0.05);
  //gStyle->SetTitleXOffset(0.1);
  gStyle->SetLabelSize(0.04,"xyz");
  gStyle->SetOptFit(0);
  //gStyle->SetStatY(0.85);
  //gStyle->SetStatX(0.975);
  //gStyle->SetStatW(.09);

  gStyle->SetFillStyle(0000); 
  //gStyle->SetStatStyle(0); 
  //gStyle->SetTitleStyle(0); 
  //gStyle->SetCanvasBorderSize(0); 
  //gStyle->SetFrameBorderSize(0); 
  //gStyle->SetLegendBorderSize(0); 
  //gStyle->SetStatBorderSize(0); 
  //gStyle->SetTitleBorderSize(0);

  ofstream params(TString::Format("%s_EQ2EtrueFitParams.dat",geom.Data()).Data()); //Output file for fit parameters
  TString fileName  = TString::Format("HistMeans_%s.dat",geom.Data());
  ifstream infile(fileName.Data());

  int type0lowOffset = 7;
  int type1lowOffset = 13;
  int type23lowOffset = 15;
  
  int type0highOffset = 2;
  int type1highOffset = 4;
  int type23highOffset = 5;

  vector <string> nameRow(23);

  vector <double> EtrueMin;
  vector <double> EtrueMax;
  vector <double> EtrueMid;
  vector <double> EQtype0E;
  vector <double> EQtype1E;
  vector <double> EQtype2E;
  vector <double> EQtype3E;
  vector <double> EQtype23E;
  vector <double> EQtype0W;
  vector <double> EQtype1W;
  vector <double> EQtype2W;
  vector <double> EQtype3W;
  vector <double> EQtype23W;
  // Errors
  vector <double> err_EQtype0E;
  vector <double> err_EQtype1E;
  vector <double> err_EQtype2E;
  vector <double> err_EQtype3E;
  vector <double> err_EQtype23E;
  vector <double> err_EQtype0W;
  vector <double> err_EQtype1W;
  vector <double> err_EQtype2W;
  vector <double> err_EQtype3W;
  vector <double> err_EQtype23W;

  double nEvents;
  int nHist,Emin,Emax;
  double type0E,type1E,type2E,type3E,type23E;
  double type0W,type1W,type2W,type3W,type23W;
  
  double err_type0E,err_type1E,err_type2E,err_type3E,err_type23E;
  double err_type0W,err_type1W,err_type2W,err_type3W,err_type23W;
  
  for (int i=0; i<23; i++) {
    infile >> nameRow[i];
    //cout << nameRow[i] << " " ;
  }
  //cout << endl;
  
  while (infile >> nHist >> Emin >> Emax >> type0E >> err_type0E >> type0W >> err_type0W >>
	 type1E >> err_type1E >> type1W >> err_type1W >> type2E >> err_type2E >> type2W >> err_type2W >>
	 type3E >> err_type3E >> type3W >> err_type3W >> type23E >> err_type23E >> type23W >> err_type23W) {
    //if (Emin<=50.) continue;
    cout << Emin << " " << type0E << endl;
    EtrueMin.push_back(Emin);
    EtrueMax.push_back(Emax);
    EtrueMid.push_back((Emin+Emax)/2.);
    EQtype0E.push_back(type0E);
    EQtype1E.push_back(type1E);
    EQtype2E.push_back(type2E);
    EQtype3E.push_back(type3E);
    EQtype23E.push_back(type23E);
    EQtype0W.push_back(type0W);
    EQtype1W.push_back(type1W);
    EQtype2W.push_back(type2W);
    EQtype3W.push_back(type3W);
    EQtype23W.push_back(type23W);

    err_EQtype0E.push_back(err_type0E);
    err_EQtype1E.push_back(err_type1E);
    err_EQtype2E.push_back(err_type2E);
    err_EQtype3E.push_back(err_type3E);
    err_EQtype23E.push_back(err_type23E);
    err_EQtype0W.push_back(err_type0W);
    err_EQtype1W.push_back(err_type1W);
    err_EQtype2W.push_back(err_type2W);
    err_EQtype3W.push_back(err_type3W);
    err_EQtype23W.push_back(err_type23W);

    if (infile.eof()) break;
  }
  infile.close();

  
  
  int numDataPoints = EtrueMin.size();
  //cout << numDataPoints << " " << EQtype3E.size() << " " <<EtrueMid.size() <<  endl;

  TF1 *func0 = new TF1("func0","[0]+[1]*x+[2]/x+[3]/(x*x)",40.,1000.);
  //TF1 *func0 = new TF1("func0","[0]+[1]*x+[2]/(x+[3])+[4]/((x+[5])*(x+[5]))",10.,1000.);
  func0->SetParameters(100.,1.,0.,0.);
  func0->SetLineColor(kBlue);

  TF1 *func1 = new TF1("func1","[0]+[1]*x+[2]/x+[3]/(x*x)",40.,1000.);
  func1->SetParameters(100.,1.,0.,0.);
  func1->SetLineColor(kRed);

  TF1 *func23 = new TF1("func23","[0]+[1]*x+[2]/x+[3]/(x*x)",40.,1000.);
  //TF1 *func23 = new TF1("func0","[0]+[1]*x+[2]/(x+[3])+[4]/((x+[5])*(x+[5]))",10.,1000.);
  func23->SetParameters(100.,1.,0.,0.);
  func23->SetLineColor(kGreen);


  TGraphErrors *t0E = new TGraphErrors(numDataPoints-(type0lowOffset+type0highOffset),&EQtype0E[type0lowOffset-1],&EtrueMid[type0lowOffset-1],&err_EQtype0E[type0lowOffset-1],0); 
  t0E->SetMarkerColor(kBlue);
  t0E->SetMarkerStyle(21);
  //t0E->SetMarkerSize(.8);
  TGraphErrors *t1E = new TGraphErrors(numDataPoints-(type1lowOffset+type1highOffset),&EQtype1E[type1lowOffset-1],&EtrueMid[type1lowOffset-1],&err_EQtype1E[type1lowOffset-1],0);
  t1E->SetMarkerColor(kRed);
  t1E->SetMarkerStyle(20);
  TGraphErrors *t2E = new TGraphErrors(numDataPoints-11,&EQtype2E[8],&EtrueMid[8],&err_EQtype2E[8],0);
  TGraphErrors *t3E = new TGraphErrors(numDataPoints-11,&EQtype3E[8],&EtrueMid[8],&err_EQtype3E[8],0);
  TGraphErrors *t23E = new TGraphErrors(numDataPoints-(type23lowOffset+type23highOffset),&EQtype23E[type23lowOffset-1],&EtrueMid[type23lowOffset-1],&err_EQtype23E[type23lowOffset-1],0);
  t23E->SetMarkerColor(kGreen);
  t23E->SetMarkerStyle(22);

  //t0E->Fit(func0,"","",EQtype0E[type0lowOffset-1],800.);
  //t1E->Fit(func1,"","",EQtype1E[type1lowOffset-1],800.);
  //t23E->Fit(func23,"","",EQtype23E[type23lowOffset-1],800.);

  t0E->Fit(func0,"R");
  t1E->Fit(func1,"R");
  t23E->Fit(func23,"R");


  //Now Calculate the extrapolations to zero
  
  //East 0
  double xE0 = EQtype0E[type0lowOffset-1];
  double alphaE0 = ( ( func0->GetParameter(1)*xE0 - func0->GetParameter(2)/xE0 - 2.*func0->GetParameter(3)/(xE0*xE0) ) / 
	      ( func0->GetParameter(0) + func0->GetParameter(1)*xE0 + func0->GetParameter(2)/xE0 + func0->GetParameter(3)/(xE0*xE0) ) );
  double cE0 =  ( func0->GetParameter(0) + func0->GetParameter(1)*xE0 + func0->GetParameter(2)/xE0 + func0->GetParameter(3)/(xE0*xE0) ) / TMath::Power(xE0,alphaE0) ;

  //East 1
  double xE1 = EQtype1E[type1lowOffset-1];
  double alphaE1 = ( ( func1->GetParameter(1)*xE1 - func1->GetParameter(2)/xE1 - 2.*func1->GetParameter(3)/(xE1*xE1) ) / 
	      ( func1->GetParameter(0) + func1->GetParameter(1)*xE1 + func1->GetParameter(2)/xE1 + func1->GetParameter(3)/(xE1*xE1) ) );
  double cE1 =  ( func1->GetParameter(0) + func1->GetParameter(1)*xE1 + func1->GetParameter(2)/xE1 + func1->GetParameter(3)/(xE1*xE1) ) / TMath::Power(xE1,alphaE1) ;

  //East 23
  double xE23 = EQtype23E[type23lowOffset-1];
  double alphaE23 = ( ( func23->GetParameter(1)*xE23 - func23->GetParameter(2)/xE23 - 2.*func23->GetParameter(3)/(xE23*xE23) ) / 
	      ( func23->GetParameter(0) + func23->GetParameter(1)*xE23 + func23->GetParameter(2)/xE23 + func23->GetParameter(3)/(xE23*xE23) ) );
  double cE23 =  ( func23->GetParameter(0) + func23->GetParameter(1)*xE23 + func23->GetParameter(2)/xE23 + func23->GetParameter(3)/(xE23*xE23) ) / TMath::Power(xE23,alphaE23) ;

  params << std::fixed << std::setprecision(17);
  
  params << "Type0E: " << func0->GetParameter(0) << " " << func0->GetParameter(1) << " " << func0->GetParameter(2) << " " << func0->GetParameter(3) 
	 << " " << xE0 << " " << alphaE0 << " " << cE0 << endl;
  params << "Type1E: " << func1->GetParameter(0) << " " << func1->GetParameter(1) << " " << func1->GetParameter(2) << " " << func1->GetParameter(3) 
	 << " " << xE1 << " " << alphaE1 << " " << cE1 << endl;
  params << "Type23E: " << func23->GetParameter(0) << " " << func23->GetParameter(1) << " " << func23->GetParameter(2) << " " << func23->GetParameter(3) 
	 << " " << xE23 << " " << alphaE23 << " " << cE23 << endl;

  TGraphErrors *t0W = new TGraphErrors(numDataPoints-(type0lowOffset+type0highOffset),&EQtype0W[type0lowOffset-1],&EtrueMid[type0lowOffset-1],&err_EQtype0W[type0lowOffset-1],0); 
  t0W->SetMarkerColor(kBlue);
  t0W->SetMarkerStyle(21);
  TGraphErrors *t1W = new TGraphErrors(numDataPoints-(type1lowOffset+type1highOffset),&EQtype1W[type1lowOffset-1],&EtrueMid[type1lowOffset-1],&err_EQtype1W[type1lowOffset-1],0);
  t1W->SetMarkerColor(kRed);
  t1W->SetMarkerStyle(20);
  TGraphErrors *t2W = new TGraphErrors(numDataPoints-11,&EQtype2W[8],&EtrueMid[8],&err_EQtype2W[8],0);
  TGraphErrors *t3W = new TGraphErrors(numDataPoints-11,&EQtype3W[8],&EtrueMid[8],&err_EQtype3W[8],0);
  TGraphErrors *t23W = new TGraphErrors(numDataPoints-(type23lowOffset+type23highOffset),&EQtype23W[type23lowOffset-1],&EtrueMid[type23lowOffset-1],&err_EQtype23W[type23lowOffset-1],0);
  t23W->SetMarkerColor(kGreen);
  t23W->SetMarkerStyle(22);

  t0W->Fit(func0,"","",EQtype0W[type0lowOffset-1],800.);
  t1W->Fit(func1,"","",EQtype1W[type1lowOffset-1],800.);
  t23W->Fit(func23,"","",EQtype23W[type23lowOffset-1],800.);

  //West 0
  double xW0 = EQtype0W[type0lowOffset-1];
  double alphaW0 = ( ( func0->GetParameter(1)*xW0 - func0->GetParameter(2)/xW0 - 2.*func0->GetParameter(3)/(xW0*xW0) ) / 
	      ( func0->GetParameter(0) + func0->GetParameter(1)*xW0 + func0->GetParameter(2)/xW0 + func0->GetParameter(3)/(xW0*xW0) ) );
  double cW0 =  ( func0->GetParameter(0) + func0->GetParameter(1)*xW0 + func0->GetParameter(2)/xW0 + func0->GetParameter(3)/(xW0*xW0) ) / TMath::Power(xW0,alphaW0) ;

  //East 1
  double xW1 = EQtype1W[type1lowOffset-1];
  double alphaW1 = ( ( func1->GetParameter(1)*xW1 - func1->GetParameter(2)/xW1 - 2.*func1->GetParameter(3)/(xW1*xW1) ) / 
	      ( func1->GetParameter(0) + func1->GetParameter(1)*xW1 + func1->GetParameter(2)/xW1 + func1->GetParameter(3)/(xW1*xW1) ) );
  double cW1 =  ( func1->GetParameter(0) + func1->GetParameter(1)*xW1 + func1->GetParameter(2)/xW1 + func1->GetParameter(3)/(xW1*xW1) ) / TMath::Power(xW1,alphaW1) ;

  //East 23
  double xW23 = EQtype23W[type23lowOffset-1];
  double alphaW23 = ( ( func23->GetParameter(1)*xW23 - func23->GetParameter(2)/xW23 - 2.*func23->GetParameter(3)/(xW23*xW23) ) / 
	      ( func23->GetParameter(0) + func23->GetParameter(1)*xW23 + func23->GetParameter(2)/xW23 + func23->GetParameter(3)/(xW23*xW23) ) );
  double cW23 =  ( func23->GetParameter(0) + func23->GetParameter(1)*xW23 + func23->GetParameter(2)/xW23 + func23->GetParameter(3)/(xW23*xW23) ) / TMath::Power(xW23,alphaW23) ;

  params << "Type0W: " << func0->GetParameter(0) << " " << func0->GetParameter(1) << " " << func0->GetParameter(2) << " " << func0->GetParameter(3) 
	 << " " << xW0 << " " << alphaW0 << " " << cW0 << endl;
  params << "Type1W: " << func1->GetParameter(0) << " " << func1->GetParameter(1) << " " << func1->GetParameter(2) << " " << func1->GetParameter(3) 
	 << " " << xW1 << " " << alphaW1 << " " << cW1 << endl;
  params << "Type23W: " << func23->GetParameter(0) << " " << func23->GetParameter(1) << " " << func23->GetParameter(2) << " " << func23->GetParameter(3) 
	 << " " << xW23 << " " << alphaW23 << " " << cW23 << endl;
 
  TMultiGraph *east = new TMultiGraph();
  east->Add(t0E);
  east->Add(t1E);
  east->Add(t23E);
  east->SetTitle("East Detector");

  TMultiGraph *west = new TMultiGraph();
  west->Add(t0W);
  west->Add(t1W);
  west->Add(t23W);
  west->SetTitle("West Detector");
  
  TCanvas *c1 = new TCanvas("c1", "c1", 1200, 500);
  c1->Divide(2,1);
  c1->cd(1);
  east->Draw("AP");
  east->SetMinimum(0.);
  east->SetMaximum(800.);
  east->GetXaxis()->SetLimits(0.,800.);
  east->GetXaxis()->SetTitle("E_{Q} (keV)");
  east->GetYaxis()->SetTitle("E_{true} (keV)");
  east->GetXaxis()->CenterTitle();
  east->GetYaxis()->CenterTitle();
  east->GetXaxis()->SetTitleOffset(1.4);
  east->GetYaxis()->SetTitleOffset(1.4);
  east->Draw("AP");

  TLegend *leg0 = new TLegend(0.55,0.25,0.85,0.375);
  leg0->AddEntry(t0E,"  Type 0","p");
  leg0->AddEntry(t1E,"  Type 1","p");
  leg0->AddEntry(t23E,"  Type 2/3","p");
  leg0->Draw();

  c1->cd(2);
  west->Draw("AP");
  west->SetMinimum(0.);
  west->SetMaximum(800.);
  west->GetXaxis()->SetLimits(0.,800.);
  west->GetXaxis()->SetTitle("E_{Q} (keV)");
  west->GetYaxis()->SetTitle("E_{true} (keV)");
  west->GetXaxis()->CenterTitle();
  west->GetYaxis()->CenterTitle();
  west->GetXaxis()->SetTitleOffset(1.4);
  west->GetYaxis()->SetTitleOffset(1.4);
  west->Draw("AP");

  leg0->Draw();

  // Draw the extrapolated functions
  c1->cd(1);

  TF1 *e0 = new TF1("e0","[0]*TMath::Power(x,[1])",0.,xE0);
  e0->SetParameters(cE0, alphaE0);
  e0->SetLineColor(kBlue);
  e0->SetLineStyle(2);
  e0->Draw("SAME");

  TF1 *e1 = new TF1("e1","[0]*TMath::Power(x,[1])",0.,xE1);
  e1->SetParameters(cE1, alphaE1);
  e1->SetLineColor(kRed);
  e1->SetLineStyle(2);
  e1->Draw("SAME");

  TF1 *e23 = new TF1("e23","[0]*TMath::Power(x,[1])",0.,xE23);
  e23->SetParameters(cE23, alphaE23);
  e23->SetLineColor(kGreen);
  e23->SetLineStyle(2);
  e23->Draw("SAME");

  c1->cd(2);

  TF1 *W0 = new TF1("W0","[0]*TMath::Power(x,[1])",0.,xW0);
  W0->SetParameters(cW0, alphaW0);
  W0->SetLineColor(kBlue);
  W0->SetLineStyle(2);
  W0->Draw("SAME");

  TF1 *W1 = new TF1("W1","[0]*TMath::Power(x,[1])",0.,xW1);
  W1->SetParameters(cW1, alphaW1);
  W1->SetLineColor(kRed);
  W1->SetLineStyle(2);
  W1->Draw("SAME");

  TF1 *W23 = new TF1("W23","[0]*TMath::Power(x,[1])",0.,xW23);
  W23->SetParameters(cW23, alphaW23);
  W23->SetLineColor(kGreen);
  W23->SetLineStyle(2);
  W23->Draw("SAME");

  params.close();
  
  

}
