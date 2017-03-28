{
  //#include <vector>
  //gStyle->SetOptFit(1111);
  gStyle->SetTitleStyle(0000);
  gStyle->SetTitleX(.25);

  ifstream infile("rawAsymmetryByOctet_0-59_AnaCh1_50mm.dat");
  char in1[10], in2[10];
  double A, err;
  vector <double> octet;
  vector <double> asym; 
  vector <double> errorY;
  vector <double> errorX;
  double Octet = 0.;

  double weightSum = 0.;
  double aveSum = 0.;

  while (infile >> Octet >> in1 >> in2) {
    A=atof(in1);
    err=atof(in2);

    cout << Octet << " " << A << " " << err << endl;

    if (TMath::IsNaN(A)) {continue;}
    
    else {
      
      octet.push_back(Octet);
      asym.push_back(A);
      errorY.push_back(err); 
      aveSum += 1./(err*err)*A;
      weightSum+=1./(err*err);
      errorX.push_back(0.);
      
      Octet = Octet+1.;
    }
  }
  double finalAsym = aveSum/weightSum;
  double totalError = 1./TMath::Sqrt(weightSum);

  for (int i=0; i<octet.size(); i++) {
    cout << octet[i] << " " << asym[i] << " " << errorY[i] << endl;
  }

  TCanvas *c1 = new TCanvas("c1");

  TGraphErrors *g = new TGraphErrors(octet.size(), &octet[0],&asym[0],&errorX[0], &errorY[0]);
  g->SetTitle("Raw Measured Asymmetry");
 
  g->SetMarkerStyle(20);
  g->SetLineWidth(2);
  g->GetXaxis()->SetLimits(-2., octet[octet.size()-1]+2.);
  

  TF1 *fit = new TF1("fit","[0]",octet[0], octet[octet.size()-1]);
  fit->SetLineColor(kRed);
  fit->SetLineWidth(3);

  g->Fit("fit","R");

  
  g->Draw("AP");
  g->SetMinimum(0.03);
  g->SetMaximum(0.07);
  c1->Update();

  cout << "final Asym = " << finalAsym << " +/- " << totalError << endl;
} 
