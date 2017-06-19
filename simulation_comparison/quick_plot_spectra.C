#include <vector>
#include <iostream>
#include <cstdlib>

vector <vector <double> > returnSourcePosition (Int_t runNumber, string src) {
  Char_t temp[500];
  sprintf(temp,"%s/source_list_%i.dat",getenv("SOURCE_LIST"),runNumber);
  ifstream file(temp);
  cout << src << endl;
  int num = 0;
  file >> num;
  cout << num << endl;
  int srcNum = 0;
  string src_check;
  for (int i=0; i<num;srcNum++,i++) {
    file >> src_check;
    cout << src_check << endl;
    if (src_check==src) break;   
  }
  cout << "The source Number is: " << srcNum << endl;
  if (srcNum==num) {
    cout << "Didn't find source in that run\n"; exit(0);
  }
  file.close();
  
  sprintf(temp,"%s/source_positions_%i.dat",getenv("SOURCE_POSITIONS"),runNumber);
  file.open(temp);
  
  vector < vector < double > > srcPos;
  srcPos.resize(2,vector <double> (3,0.));
  
  for (int i=0; i<srcNum+1; i++) {
    for (int j=0; j<2; j++) {
      for (int jj=0; jj<3; jj++) {
	file >> srcPos[j][jj];
      }
    }
  }
  return srcPos;
}

void quick_plot_spectra(Int_t runNumber, string src) {
  Char_t temp[500];
  gStyle->SetOptStat(0);
  sprintf(temp,"%s/replay_pass3_%i.root",getenv("REPLAY_PASS3"),runNumber);
  TFile *data = new TFile(temp,"READ");
  TTree *TData = (TTree*)(data->Get("pass3"));
  sprintf(temp,"%s/sources/revCalSim_%i_%s.root",getenv("REVCALSIM"),runNumber, src.c_str());
  TFile *sim = new TFile(temp,"READ");
  TTree *TSim = (TTree*)(sim->Get("revCalSim"));
  cout << "Opened Files\n";

  string srcShort = src;
  srcShort.erase(2);
  vector < vector <double> > srcPos = returnSourcePosition(runNumber,srcShort);
  
  for (int j=0; j<2; j++) {
      for (int jj=0; jj<3; jj++) {
	cout << srcPos[j][jj] << " ";
      }
      cout << endl;
    }
  
  cout << "Got source position\n";
  sprintf(temp,"%s East",src.c_str());
  TH1D *data_E_0 = new TH1D("data_E_0", temp, 400, 0., 1200.);
  TH1D *sim_E_0 = new TH1D("sim_E_0", "East 0 Evts", 400, 0., 1200.);
  sim_E_0->SetLineColor(kRed);
  sprintf(temp,"%s West",src.c_str());
  TH1D *data_W_0 = new TH1D("data_W_0", temp, 400, 0., 1200.);
  TH1D *sim_W_0 = new TH1D("sim_W_0", "West 0 Evts", 400, 0., 1200.);
  sim_W_0->SetLineColor(kRed);
  Char_t cuts[1000];
  
  //0 types

  Double_t integralE, integralW, scaleValE, scaleValW;

  Double_t normMin = 0., normMax = 1200.;

  normMin = src=="Bi207" ? 900. : src=="Sn113" ? 300. : 80.;
  normMax = src=="Bi207" ? 1100. : src=="Sn113" ? 450. : 170.;

  TCanvas *c1 = new TCanvas("c1","c1", 1200, 400);
  c1->Divide(2,1);
  c1->cd(1);
  sprintf(cuts,"PID==1 && Type==0 && Side==0 && EvisE>0. && xE.center>(%f-2.*%f) && xE.center<(%f+2.*%f) && yE.center>(%f-2.*%f) && yE.center<(%f+2.*%f)",srcPos[0][0],srcPos[0][2],srcPos[0][0],srcPos[0][2],srcPos[0][1],srcPos[0][2],srcPos[0][1],srcPos[0][2]);  
  TData->Draw("Erecon>>data_E_0", cuts);

  data_E_0->GetXaxis()->SetRangeUser(normMin,normMax);
  integralE = data_E_0->Integral();
  //Double_t
  data_E_0->GetXaxis()->SetRangeUser(0.,1200.);

  sprintf(cuts,"PID==1 && side==0 && type==0 && EvisE>0.");
  TSim->Draw("Erecon>>sim_E_0",cuts,"SAME");
  sim_E_0->GetXaxis()->SetRangeUser(normMin,normMax);
  scaleValE = integralE/sim_E_0->Integral();
  //Double_t scaleValE = data_E_0->GetBinContent(data_E_0->GetMaximumBin())/sim_E_0->GetBinContent(sim_E_0->GetMaximumBin());
  sim_E_0->GetXaxis()->SetRangeUser(0.,1200.);
  sim_E_0->Scale(scaleValE);
  sim_E_0->Draw("SAME");

  //TCanvas *c2 = new TCanvas("c2");
  c1->cd(2);

  sprintf(cuts,"PID==1 && Type==0 && Side==1 && EvisW>0. && xW.center>(%f-2.*%f) && xW.center<(%f+2.*%f) && yW.center>(%f-2.*%f) && yW.center<(%f+2.*%f)",srcPos[1][0],srcPos[1][2],srcPos[1][0],srcPos[1][2],srcPos[1][1],srcPos[1][2],srcPos[1][1],srcPos[1][2]);  
  TData->Draw("Erecon>>data_W_0", cuts);
  data_W_0->GetXaxis()->SetRangeUser(normMin,normMax);
  integralW = data_W_0->Integral();
  data_W_0->GetXaxis()->SetRangeUser(0.,1200.);

  sprintf(cuts,"PID==1 && side==1 && type==0 && EvisW>0.");
  TSim->Draw("Erecon>>sim_W_0",cuts,"SAME");
  sim_W_0->GetXaxis()->SetRangeUser(normMin,normMax);
  scaleValW = integralW/sim_W_0->Integral();
  //Double_t scaleValW = data_W_0->GetBinContent(data_W_0->GetMaximumBin())/sim_W_0->GetBinContent(sim_W_0->GetMaximumBin());

  sim_W_0->GetXaxis()->SetRangeUser(0.,1200.);
  sim_W_0->Scale(scaleValW);
  sim_W_0->Draw("SAME");

  
}
