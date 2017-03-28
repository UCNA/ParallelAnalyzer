

void find23separation(TString geom) {

  gStyle->SetOptStat(0);
  //gStyle->SetStatFontSize(0.030);
  //gStyle->SetOptFit(1111);
  //gStyle->SetOptTitle(0);
  //gStyle->SetTitleSize(0.08,"t");
  //gStyle->SetTitleX(0.17);
  //gStyle->SetTitleAlign(13);
  gStyle->SetTitleOffset(1.30, "x");
  gStyle->SetTitleOffset(1.1, "y");
  //gStyle->SetPadTickX(1);
  //gStyle->SetPadTickY(1);
  //gStyle->SetNdivisions(510,"X");
  //gStyle->SetNdivisions(510,"Y");
  //gStyle->SetNdivisions(9,"Z");
  gStyle->SetPadLeftMargin(0.13); // 0.13
  //gStyle->SetPadRightMargin(0.15); // 0.04
  gStyle->SetPadBottomMargin(0.17); // 0.30
  gStyle->SetLabelSize(0.045, "XYZ");
  //gStyle->SetLabelSize(0.045, "Y");
  //gStyle->SetLabelSize(0.045, "Z");
  //gStyle->SetLabelOffset(0.00, "X");
  //gStyle->SetLabelOffset(0.01, "Y");
  gStyle->SetTitleSize(0.050, "X");
  gStyle->SetTitleSize(0.050, "Y");

  //Borders on Legends and titles and stats
  gStyle->SetFillStyle(0000); 
  //gStyle->SetStatStyle(0); 
  //gStyle->SetTitleStyle(0); 
  //gStyle->SetCanvasBorderSize(0); 
  //gStyle->SetFrameBorderSize(0); 
  gStyle->SetLegendBorderSize(0); 
  //gStyle->SetStatBorderSize(0); 
  //gStyle->SetTitleBorderSize(0);
  
  int octMin = 0;
  int octMax = 0;
  double sepVal = 6.6;

  static const  Int_t numEnergyBins = 8;
  Double_t energyStart = 0.;
  Double_t energyBinWidth = 100.;
  
  TH1D *hType2E[numEnergyBins];
  TH1D *hType3E[numEnergyBins];
  TH1D *hType2W[numEnergyBins];
  TH1D *hType3W[numEnergyBins];

  TH1D *hScint2E[numEnergyBins];
  TH1D *hScint3E[numEnergyBins];
  TH1D *hScint2W[numEnergyBins];
  TH1D *hScint3W[numEnergyBins];
  
  if ( geom==TString("2011-2012") ) octMin = 0, octMax = 59;
  else if ( geom==TString("2012-2013") ) octMin = 80, octMax = 121;
  else if ( geom==TString("2012-2013_isobutane") ) octMin = 67, octMax = 79, sepVal = 5.4 ; 
  else { std::cout << "Bad geometry\n"; exit(0); }

  //TH1D::AddDirectory(false);
  
  TFile *f = new TFile(TString::Format("Backscatters_%i-%i.root",octMin,octMax),
		       "UPDATE");

  for ( Int_t hist=0; hist<numEnergyBins; ++hist ) {

    Double_t binLowEdge = hist*energyBinWidth + energyStart;
    Double_t binHighEdge = binLowEdge + energyBinWidth;
    
    hType2E[hist] = (TH1D*)f->Get(TString::Format("hType2E_%0.0f-%0.0f",binLowEdge,binHighEdge));//->Clone();
    hType3E[hist] = (TH1D*)f->Get(TString::Format("hType3E_%0.0f-%0.0f",binLowEdge,binHighEdge));//->Clone();
    hType2W[hist] = (TH1D*)f->Get(TString::Format("hType2W_%0.0f-%0.0f",binLowEdge,binHighEdge));//->Clone();
    hType3W[hist] = (TH1D*)f->Get(TString::Format("hType3W_%0.0f-%0.0f",binLowEdge,binHighEdge));//->Clone();
    
    hScint2E[hist] = (TH1D*)f->Get(TString::Format("hScint2E_%0.0f-%0.0f",binLowEdge,binHighEdge));//->Clone();
    hScint3E[hist] = (TH1D*)f->Get(TString::Format("hScint3E_%0.0f-%0.0f",binLowEdge,binHighEdge));//->Clone();
    hScint2W[hist] = (TH1D*)f->Get(TString::Format("hScint2W_%0.0f-%0.0f",binLowEdge,binHighEdge));//->Clone();
    hScint3W[hist] = (TH1D*)f->Get(TString::Format("hScint3W_%0.0f-%0.0f",binLowEdge,binHighEdge));//->Clone();
  }
  
  //delete f;

  //Now I need to make new histograms of misidentified fractions
  // and fit the minimum region to find the minimum point

  Int_t nbins = 100;
  Double_t binWidth = 20. / nbins;
  TH1D *misIDE[numEnergyBins];
  Int_t total2E[numEnergyBins];
  Int_t total3E[numEnergyBins];
  Int_t totalE[numEnergyBins];
  TF1 *funcE[numEnergyBins];

  Double_t enBinMidE[numEnergyBins], fitValE[numEnergyBins], fitValErrE[numEnergyBins], maxBinValE[numEnergyBins];


  for ( Int_t hist=0; hist<numEnergyBins; ++hist ) {


    Double_t binLowEdge = hist*energyBinWidth + energyStart;
    Double_t binHighEdge = binLowEdge + energyBinWidth;  
    
    misIDE[hist] = new TH1D(TString::Format("misIDE_%0.0f-%0.0f",binLowEdge,binHighEdge),
			   TString::Format("Properly Identified Fraction %0.0f-%0.0f keV",binLowEdge,binHighEdge),
			   nbins, binWidth/2., 20.+binWidth/2.);
    misIDE[hist]->GetXaxis()->SetTitle("Position of Separation Cut in E_{MWPC} (keV)");
    misIDE[hist]->GetYaxis()->SetTitle("Fraction of Correctly Identified Events");
    misIDE[hist]->GetXaxis()->CenterTitle();
    misIDE[hist]->GetYaxis()->CenterTitle();
	
    total2E[hist] = hType2E[hist]->Integral(1,nbins);
    total3E[hist] = hType3E[hist]->Integral(1,nbins);
    totalE[hist] = total2E[hist] + total3E[hist];
    
    for ( Int_t cut=1; cut<=nbins; ++cut ) {
      
      Int_t corr2 = hType2E[hist]->Integral(1,cut);
      Int_t corr3 = hType3E[hist]->Integral(cut,nbins);
      Double_t frac = (Double_t)(corr2 + corr3)/(Double_t)totalE[hist];
      misIDE[hist]->SetBinContent( cut, frac );
    }
    
    Double_t max = misIDE[hist]->GetXaxis()->GetBinCenter( misIDE[hist]->GetMaximumBin() );
    
    funcE[hist] = new TF1(TString::Format("funcE%i",hist),"[0]+[1]*(x-[3])+[2]*(x-[3])**2",max-1.,max+1.);
    funcE[hist]->SetParameters(1.,1.,-1.,max);
    //func[hist] = new TF1(TString::Format("func%i",hist),"[0]*TMath::Gaus(x,[1],[2])",max-1.,max+1.);
    //funcE[hist]->SetParameters(1.,1.,1.);
    
    misIDE[hist]->Fit(funcE[hist],"R");
    
    /*Double_t maxFit = ( ( -funcE[hist]->GetParameter(2) +
			TMath::Sqrt( TMath::Power(funcE[hist]->GetParameter(2),2)-3.*funcE[hist]->GetParameter(1)*funcE[hist]->GetParameter(3)) )
			/ funcE[hist]->GetParameter(1) ) + funcE[hist]->GetParameter(4);//

    if ( maxFit < (max-1.) || maxFit > (max+1.) ) {
      maxFit = ( ( -funcE[hist]->GetParameter(2) -
		   TMath::Sqrt( TMath::Power(funcE[hist]->GetParameter(2),2)-3.*funcE[hist]->GetParameter(1)*funcE[hist]->GetParameter(3)) )
		 / funcE[hist]->GetParameter(1) )  + funcE[hist]->GetParameter(4);
		 }*/

    Double_t maxFit = -funcE[hist]->GetParameter(1) / (2.*funcE[hist]->GetParameter(2)) + funcE[hist]->GetParameter(3);
    
    std::cout << "Max Efficiency Bin Center = " << max << std::endl;
    std::cout << "Max Efficiency Fit Val = " << maxFit << std::endl;
    
    enBinMidE[hist] = binHighEdge - energyBinWidth/2.;
    fitValE[hist] = maxFit;
    fitValErrE[hist] = TMath::Sqrt( TMath::Power( funcE[hist]->GetParError(1) / (2.*funcE[hist]->GetParameter(2)) , 2 ) +
				   TMath::Power( funcE[hist]->GetParError(2) *
    						 ( funcE[hist]->GetParameter(1) / (2.*TMath::Power(funcE[hist]->GetParameter(2),2) ) ), 2 ) );
    maxBinValE[hist] = max;
  }

  TH1D *misIDW[numEnergyBins];
  Int_t total2W[numEnergyBins];
  Int_t total3W[numEnergyBins];
  Int_t totalW[numEnergyBins];
  TF1 *funcW[numEnergyBins];

  Double_t enBinMidW[numEnergyBins], fitValW[numEnergyBins], fitValErrW[numEnergyBins], maxBinValW[numEnergyBins];


  for ( Int_t hist=0; hist<numEnergyBins; ++hist ) {


    Double_t binLowEdge = hist*energyBinWidth + energyStart;
    Double_t binHighEdge = binLowEdge + energyBinWidth;  
    
    misIDW[hist] = new TH1D(TString::Format("misIDW_%0.0f-%0.0f",binLowEdge,binHighEdge),
			   TString::Format("Properly Identified Fraction %0.0f-%0.0f keV",binLowEdge,binHighEdge),
			   nbins, binWidth/2., 20.+binWidth/2.);
    misIDW[hist]->GetXaxis()->SetTitle("Position of Separation Cut in E_{MWPC} (keV)");
    misIDW[hist]->GetYaxis()->SetTitle("Fraction of Correctly Identified Events");
    misIDW[hist]->GetXaxis()->CenterTitle();
    misIDW[hist]->GetYaxis()->CenterTitle();
	
    total2W[hist] = hType2W[hist]->Integral(1,nbins);
    total3W[hist] = hType3W[hist]->Integral(1,nbins);
    totalW[hist] = total2W[hist] + total3W[hist];
    
    for ( Int_t cut=1; cut<=nbins; ++cut ) {
      
      Int_t corr2 = hType2W[hist]->Integral(1,cut);
      Int_t corr3 = hType3W[hist]->Integral(cut,nbins);
      Double_t frac = (Double_t)(corr2 + corr3)/(Double_t)totalW[hist];
      misIDW[hist]->SetBinContent( cut, frac );
    }
    
    Double_t max = misIDW[hist]->GetXaxis()->GetBinCenter( misIDW[hist]->GetMaximumBin() );
    
    funcW[hist] = new TF1(TString::Format("funcW%i",hist),"[0]+[1]*(x-[3])+[2]*(x-[3])**2",max-1.,max+1.);
    funcW[hist]->SetParameters(1.,1.,-1,max);
    //func[hist] = new TF1(TString::Format("func%i",hist),"[0]*TMath::Gaus(x,[1],[2])",max-1.,max+1.);
    //func[hist]->SetParameters(1.,1.,1.);
    
    misIDW[hist]->Fit(funcW[hist],"R");

    /*Double_t maxFit = ( ( -funcW[hist]->GetParameter(2) +
			TMath::Sqrt( TMath::Power(funcW[hist]->GetParameter(2),2)-3.*funcW[hist]->GetParameter(1)*funcW[hist]->GetParameter(3)) )
			/ funcW[hist]->GetParameter(1) ) + funcW[hist]->GetParameter(4);//-funcW[hist]->GetParameter(1) / (2.*funcW[hist]->GetParameter(2));

    if ( maxFit < (max-1.) || maxFit > (max+1.) ) {
      maxFit = ( ( -funcW[hist]->GetParameter(2) -
		   TMath::Sqrt( TMath::Power(funcW[hist]->GetParameter(2),2)-3.*funcW[hist]->GetParameter(1)*funcW[hist]->GetParameter(3)) )
		 / funcW[hist]->GetParameter(1) )  + funcW[hist]->GetParameter(4);
		 }*/

    Double_t maxFit = -funcW[hist]->GetParameter(1) / (2.*funcW[hist]->GetParameter(2)) + funcW[hist]->GetParameter(3);

    std::cout << "Max Efficiency Bin Center = " << max << std::endl;
    std::cout << "Max Efficiency Fit Val = " << maxFit << std::endl;
    
    enBinMidW[hist] = binHighEdge - energyBinWidth/2.;
    fitValW[hist] = maxFit;
    fitValErrW[hist] = TMath::Sqrt( TMath::Power( funcW[hist]->GetParError(1) / (2.*funcW[hist]->GetParameter(2)) , 2 ) +
				   TMath::Power( funcW[hist]->GetParError(2) *
    						 ( funcW[hist]->GetParameter(1) / (2.*TMath::Power(funcW[hist]->GetParameter(2),2) ) ), 2 ) );
    maxBinValW[hist] = max;
  }

  TCanvas *c1 = new TCanvas("c1","c1");
  misIDE[3]->Draw();

  TCanvas *c2 = new TCanvas("c2");


  TGraph *gE = new TGraph(numEnergyBins,enBinMidE,fitValE);
  gE->SetMarkerStyle(26);
  gE->SetTitle("Type 2/3 Separation Parameterization");
  gE->GetXaxis()->SetTitle("Scintillator Energy (keV)");
  gE->GetYaxis()->SetTitle("Wirechamber Energy Cut (keV)");
  gE->SetLineWidth(2);
  gE->GetXaxis()->CenterTitle();
  gE->GetYaxis()->CenterTitle();
  gE->SetMarkerColor(kBlue);
  gE->SetMarkerSize(1.5);
  gE->SetMaximum(12.);
  gE->Draw("AP");

  TGraph *gW = new TGraph(numEnergyBins,enBinMidW,fitValW);
  gW->SetMarkerStyle(24);
  gW->SetTitle("Type 2/3 Separation Parameterization");
  gW->GetXaxis()->SetTitle("Scintillator Energy (keV)");
  gW->GetYaxis()->SetTitle("Wirechamber Energy Cut (keV)");
  gW->GetXaxis()->CenterTitle();
  gW->GetYaxis()->CenterTitle();
  gW->SetMarkerColor(kRed);
  gW->SetMarkerSize(1.5);
  gW->Draw("PSAME");

  Double_t aveFit[numEnergyBins];
  for ( int i=0; i<numEnergyBins; ++i ) {
    aveFit[i] = ( fitValW[i] + fitValE[i] ) / 2.;
  }

  TGraph *ave = new TGraph(numEnergyBins,enBinMidW,aveFit);
  ave->SetMarkerStyle(1);
  ave->SetTitle("Type 2/3 Separation Parameterization");
  ave->GetXaxis()->SetTitle("Scintillator Energy (keV)");
  ave->GetYaxis()->SetTitle("Wirechamber Energy Cut (keV)");
  ave->GetXaxis()->CenterTitle();
  ave->GetYaxis()->CenterTitle();
  ave->SetMarkerColor(0);
  ave->SetMarkerSize(0);
  ave->Draw("PSAME");

  Double_t fitMin = geom==TString("2012-2013_isobutane")?55.:50.;
  
  TF1 *func3 = new TF1("func3","[0]+[1]*TMath::Exp(-x/[2])",fitMin, 900.);
  func3->SetParameters(4,5.4,200);
  func3->SetLineColor(1);
  func3->SetLineStyle(7);
  func3->SetLineWidth(2);

  ave->Fit(func3,"R");

  TLegend *legFit = new TLegend(0.5,0.6,0.8,0.8);
  legFit->AddEntry(gE,"East","p");
  legFit->AddEntry(gW,"West","p");
  legFit->AddEntry(func3,"Fit to average");
  legFit->Draw("SAME");
  
  c2->Print(TString::Format("EreconVsCut_%s.pdf",geom.Data()));

  std::ofstream ofile(TString::Format("%s/backscSepParameters_%s.dat",
  				      getenv("MWPC_CALIBRATION"),geom.Data()).Data());
  //std::ofstream ofile(TString::Format("backscSepParameters_%s.dat",geom.Data()).Data());
  ofile << std::setprecision(7);
  ofile << "#p0\t\tp1\t\tp2\n"
	<< func3->GetParameter(0) << "\t" << func3->GetParameter(1) << "\t" << func3->GetParameter(2);

  ofile.close();

  f->Write();

  // Now make a TGraphErrors

  //Now we plot all
  
  ////////////////////////////////////////////////////////
  
   TCanvas *c0;
   TLegend *leg;
   TLine *line;

  for ( int i=0; i<numEnergyBins; ++i ) {
    Double_t binLowEdge = i*energyBinWidth + energyStart;
    Double_t binHighEdge = binLowEdge + energyBinWidth;  
    
    
    c0 = new TCanvas("c0","c0",1200,800);
    c0->Divide(2,1);

    c0->cd(1);

    hType2E[i]->SetTitle(TString::Format("East E_{MWPC} (%0.0f-%0.0f keV bin)",binLowEdge,binHighEdge));
    hType2E[i]->GetXaxis()->SetTitle("E_{MWPC} (keV)");
    hType2E[i]->GetXaxis()->CenterTitle();
    hType2E[i]->SetLineColor(kBlue);
    hType2E[i]->SetLineWidth(2);
    hType2E[i]->Draw();
    
    hType3E[i]->SetLineColor(kBlue);
    hType3E[i]->SetLineWidth(2);
    hType3E[i]->SetLineStyle(7);
    hType3E[i]->Draw("SAME");
    c0->Update();
    
    leg = new TLegend(0.55,0.70,0.85,0.85);
    leg->AddEntry(hType2E[i],"Type 2");
    leg->AddEntry(hType3E[i],"Type 3");
    leg->Draw("SAME");
    
    line = new TLine(fitValE[i],0.,fitValE[i],gPad->GetUymax());
    line->SetLineWidth(2);
    line->SetLineStyle(7);
    line->Draw("SAME");
    
    c0->cd(2);
    
    
    hType2W[i]->SetTitle(TString::Format("West E_{MWPC} (%0.0f-%0.0f keV)",binLowEdge,binHighEdge));
    hType2W[i]->GetXaxis()->SetTitle("E_{MWPC} (keV)");
    hType2W[i]->GetXaxis()->CenterTitle();
    hType2W[i]->SetLineColor(kBlue);
    hType2W[i]->SetLineWidth(2);
    hType2W[i]->Draw();
    
    hType3W[i]->SetLineColor(kBlue);
    hType3W[i]->SetLineWidth(2);
    hType3W[i]->SetLineStyle(7);
    hType3W[i]->Draw("SAME");
    c0->Update();
    
    leg = new TLegend(0.55,0.70,0.85,0.85);
    leg->AddEntry(hType2W[i],"Type 2");
    leg->AddEntry(hType3W[i],"Type 3");
    leg->Draw("SAME");
    
    line = new TLine(fitValW[i],0.,fitValW[i],gPad->GetUymax());
    line->SetLineWidth(2);
    line->SetLineStyle(7);
    line->Draw("SAME");

    c0->Print(TString::Format("sepPlots_%s.pdf%s",geom.Data(),i==0?"(":( i==(numEnergyBins-1)? ")" :"")));

    delete c0;
    delete leg;
    delete line;
  }

  
  
  

}
