
void plotSuperSums(int octetStart, int octetEnd) 
{
  gStyle->SetOptStat(0);
  

  TString normType = "ALL";

  Double_t normLow = 200.;
  Double_t normHigh = 700.;
  
  Double_t xAxisMax = 1200.;


  TString fileBase = "superSumPlots/SuperSum_octets_";
  fileBase+=octetStart;
  fileBase+="-";
  fileBase+=octetEnd;
  fileBase+="_AsymOff";
  
  TString fileName = TString::Format("%s.root",fileBase.Data());
  TFile *f = new TFile(fileName,"READ");

  
  TH1D *ukALL = (TH1D*)f->Get("EreconALL_uk");
  TH1D *simALL = (TH1D*)f->Get("EreconALL_sim");
  TH1D *uk0 = (TH1D*)f->Get("Erecon0_uk");
  TH1D *sim0 = (TH1D*)f->Get("Erecon0_sim");
  TH1D *uk1 = (TH1D*)f->Get("Erecon1_uk");
  TH1D *sim1 = (TH1D*)f->Get("Erecon1_sim");
  TH1D *uk23 = (TH1D*)f->Get("Erecon23_uk");
  TH1D *sim23 = (TH1D*)f->Get("Erecon23_sim");

  //Normalize
  Double_t normFactor = 1.;
  
  if (normType==TString("ALL"))
  normFactor = ukALL->Integral(ukALL->GetXaxis()->FindFixBin(normLow),ukALL->GetXaxis()->FindFixBin(normHigh))/simALL->Integral(simALL->GetXaxis()->FindFixBin(normLow),simALL->GetXaxis()->FindFixBin(normHigh));

  if (normType==TString("0"))
  normFactor = uk0->Integral(uk0->GetXaxis()->FindFixBin(normLow),uk0->GetXaxis()->FindFixBin(normHigh))/sim0->Integral(sim0->GetXaxis()->FindFixBin(normLow),sim0->GetXaxis()->FindFixBin(normHigh));

  if (normType==TString("1"))
  normFactor = uk1->Integral(uk1->GetXaxis()->FindFixBin(normLow),uk1->GetXaxis()->FindFixBin(normHigh))/sim1->Integral(sim1->GetXaxis()->FindFixBin(normLow),sim1->GetXaxis()->FindFixBin(normHigh));

  if (normType==TString("23"))
  normFactor = uk23->Integral(uk23->GetXaxis()->FindFixBin(normLow),uk23->GetXaxis()->FindFixBin(normHigh))/sim23->Integral(sim23->GetXaxis()->FindFixBin(normLow),sim23->GetXaxis()->FindFixBin(normHigh));
  
  simALL->Scale(normFactor);

  //Residuals...
  TCanvas *c1 = new TCanvas("c1","c1", 1600., 600.);
  c1->Divide(2,1);
  c1->cd(2);
  Int_t nBins = ukALL->GetNbinsX();
  Double_t min = ukALL->GetXaxis()->GetBinLowEdge(ukALL->GetXaxis()->GetFirst());
  Double_t max = ukALL->GetXaxis()->GetBinUpEdge(ukALL->GetXaxis()->GetLast());

  //residual
  TH1D *resid = new TH1D("resid","Residuals: MC-Data", nBins, min, max);
  resid->Add(simALL,ukALL,100,-100);
  resid->GetXaxis()->SetRangeUser(0., xAxisMax);
  resid->GetYaxis()->SetTitle("event rate (mHz/keV)");
  resid->SetMaximum(1.8);
  resid->SetMinimum(-1.8);
  resid->SetLineWidth(2);

  TH1D *perc_resid = new TH1D("perc_resid","Fractional Residuals: (MC-Data)/MC", nBins, min, max);
  perc_resid->Add(simALL,ukALL,100,-100);
  perc_resid->Divide(perc_resid,simALL,1.,100.);
  perc_resid->GetXaxis()->SetRangeUser(0., xAxisMax);
  //perc_resid->GetYaxis()->SetTitle("event rate (mHz/keV)");
  perc_resid->SetMaximum(.1);
  perc_resid->SetMinimum(-.1);
  perc_resid->SetLineWidth(2);


  resid->Draw();
  c1->Update();
  TLine *l = new TLine(min, 0., xAxisMax, 0.);
  l->SetLineStyle(8);
  l->Draw();

  c1->cd(1);
  
  ukALL->SetMarkerColor(kBlue);
  ukALL->SetLineColor(kBlue);
  //ukALL->SetFillStyle(3002);
  //ukALL->SetFillColor(kBlue);
  ukALL->SetLineWidth(3);
  simALL->SetMarkerColor(kRed);
  simALL->SetLineColor(kRed);
  simALL->SetMarkerSize(0.75);
  simALL->SetMarkerStyle(20);
  ukALL->GetXaxis()->SetRangeUser(0., xAxisMax);
  ukALL->GetYaxis()->SetTitle("event rate (mHz/keV)");
  ukALL->Scale(100.);
  simALL->Scale(100.);
  ukALL->SetMaximum(ukALL->GetMaximum()*1.2);
  simALL->GetXaxis()->SetRangeUser(0., xAxisMax);
  ukALL->Draw("HISTE0");
  simALL->Draw("SAMEE0");


  TCanvas *c2 = new TCanvas("c2", "c2", 1600., 600.);
  c2->Divide(3,1);
  c2->cd(1);
  
  

  uk0->SetMarkerColor(kBlue);
  uk0->SetLineColor(kBlue);
  //uk0->SetFillStyle(3002);
  //uk0->SetFillColor(kBlue);
  uk0->SetLineWidth(3);
  sim0->SetMarkerColor(kRed);
  sim0->SetLineColor(kRed);
  sim0->SetMarkerSize(0.75);
  sim0->SetMarkerStyle(20);
  uk0->GetXaxis()->SetRangeUser(0., xAxisMax);
  uk0->GetYaxis()->SetTitle("event rate (mHz/keV)");
  uk0->Scale(100.);
  sim0->Scale(100.*normFactor);
  uk0->SetMaximum(uk0->GetMaximum()*1.2);
  sim0->GetXaxis()->SetRangeUser(0., xAxisMax);
  uk0->Draw("HISTE0");
  sim0->Draw("SAMEE0");
  
  c2->cd(2);
  

  uk1->SetMarkerColor(kBlue);
  uk1->SetLineColor(kBlue);
  //uk1->SetFillStyle(3002);
  //uk1->SetFillColor(kBlue);
  uk1->SetLineWidth(3);
  uk1->SetMinimum(0.);
  sim1->SetMarkerColor(kRed);
  sim1->SetLineColor(kRed);
  sim1->SetMarkerSize(0.75);
  sim1->SetMarkerStyle(20);
  uk1->GetXaxis()->SetRangeUser(0., xAxisMax);
  uk1->GetYaxis()->SetTitle("event rate (mHz/keV)");
  uk1->Scale(100.);
  sim1->Scale(100.*normFactor);
  uk1->SetMaximum(uk1->GetMaximum()*1.2);
  sim1->GetXaxis()->SetRangeUser(0., xAxisMax);
  uk1->Draw("HISTE0");
  sim1->Draw("SAMEE0");
  

  c2->cd(3);

  TH1D *uk23 = (TH1D*)f->Get("Erecon23_uk");
  TH1D *sim23 = (TH1D*)f->Get("Erecon23_sim");
  
  uk23->SetMarkerColor(kBlue);
  uk23->SetLineColor(kBlue);
  //uk23->SetFillStyle(3002);
  //uk23->SetFillColor(kBlue);
  uk23->SetLineWidth(3);
  uk23->SetMinimum(0.);
  sim23->SetMarkerColor(kRed);
  sim23->SetLineColor(kRed);
  sim23->SetMarkerSize(0.75);
  sim23->SetMarkerStyle(20);
  uk23->GetXaxis()->SetRangeUser(0., xAxisMax);
  uk23->GetYaxis()->SetTitle("event rate (mHz/keV)");
  uk23->Scale(100.);
  sim23->Scale(100.*normFactor);
  uk23->SetMaximum(uk23->GetMaximum()*1.2);
  sim23->GetXaxis()->SetRangeUser(0., xAxisMax);
  uk23->Draw("HISTE0");
  sim23->Draw("SAMEE0");


  TString pdfFile = TString::Format("%s_normType%s_%0.1f-%0.1f.pdf",fileBase.Data(),normType.Data(), normLow, normHigh);
  //TString pdfFileStart = pdfFile + TString("(");
  //TString pdfFileEnd = pdfFile + TString
  c1->Print(TString::Format("%s[",pdfFile.Data()));
  c1->Print(pdfFile);
  c2->Print(pdfFile);
  c2->Print(TString::Format("%s]",pdfFile.Data()));
  
}
