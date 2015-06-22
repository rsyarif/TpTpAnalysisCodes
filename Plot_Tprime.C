{
  TFile *fs = new TFile("rootfiles/TprimeTToTH_M-1000_pythia8_R08_r020_KtMj_HiggsWin30_tagr07_fake02_BtagALLMj_MjIVF020_mc_subjets.root");
  TDirectory *ds = (TDirectory*) fs->GetDirectory("btaganaFatJets");
  TTree *ts = (TTree*) ds->Get("ttree");

  TFile *fb = new TFile("rootfiles/ZPrimeToTTJets_M1000GeV_R08_r020_KtMj_HiggsWin30_tagr07_fake02_BtagALLMj_MjIVF020_mc_subjets.root");
  TDirectory *db = (TDirectory*) fb->GetDirectory("btaganaFatJets");
  TTree *tb = (TTree*) db->Get("ttree");

  gStyle->SetOptStat("emrou");

  //chi wit cut

  TCanvas *c1 = new TCanvas("c1","c1",600,600);

  TH1D *h1 = new TH1D("h1","h1",40,-26,-6);
  TH1D *h2 = new TH1D("h2","h2",40,-26,-6);

  ts->Draw("log(FatJetInfo.Jet_SD_chi)>>h1","FatJetInfo.Jet_SD_chi>0&&FatJetInfo.Jet_SD_nBtagMicrojets>1");
  tb->Draw("log(FatJetInfo.Jet_SD_chi)>>h2","FatJetInfo.Jet_SD_chi>0&&FatJetInfo.Jet_SD_nBtagMicrojets>1");

  h1->SetLineColor(kBlue);
  h2->SetLineColor(kRed);

  h1->SetTitle("TprimetH M1000");
  h2->SetTitle("Zprimett M1000");

  h1->Scale(1/h1->GetEntries());
  h2->Scale(1/h2->GetEntries());

  THStack *h = new THStack("h","");
  h->Add(h1,"sames");
  h->Add(h2,"sames");
  h->Draw("nostack");
  h->GetXaxis()->SetTitle("Log(#chi)");
  h->SetTitle("cut: # of b-tagged Microjets > 1");

  TLegend *leg = new TLegend(0.2,0.65,0.4,0.85);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(h1,"T'->tH M1000","L");
  leg->AddEntry(h2,"Z'->ttbar ","L");
  leg->Draw("SAME");

  gPad->Update();

  TPaveStats *tps1 = (TPaveStats*) h1->FindObject("stats");
  tps1->SetTextColor(kBlue);
  double X1 = tps1->GetX1NDC();
  double Y1 = tps1->GetY1NDC();
  double X2 = tps1->GetX2NDC();
  double Y2 = tps1->GetY2NDC();
  TPaveStats *tps2 = (TPaveStats*) h2->FindObject("stats");
  tps2->SetTextColor(kRed);
  tps2->SetX1NDC(X1);
  tps2->SetX2NDC(X2);
  tps2->SetY1NDC(Y1-(Y2-Y1));
  tps2->SetY2NDC(Y1);
  X1 = tps2->GetX1NDC();
  Y1 = tps2->GetY1NDC();
  X2 = tps2->GetX2NDC();
  Y2 = tps2->GetY2NDC();

  gPad->Update();

  c1->SaveAs("plots/c1.pdf");

  //ROC

  const int xbin = 40; double xmin = -26; double xmax = -6;

  double denom1 = 10174; //signal, with chi>0 cut
  double denom2 = 1244; //bkg, with chi>0 cut

  double num1[xbin];
  double num2[xbin];

  double eff1[xbin];
  double eff2[xbin];

  double x[xbin];
  double inc = (xmax-xmin)/xbin;

  for(int i=0; i<xbin;i++){
    num1[i]=h1->Integral(i+1,xbin+1);
    num2[i]=h2->Integral(i+1,xbin+1);
    eff1[i]=num1[i];
    eff2[i]=num2[i];
    x[i] = xmin + inc*i;
    cout<<i+1<<". x = "<<x[i]<< ", h1 =  "<< num1[i]<< ", h2 = "<< num2[i] << endl;
  }

  double rej2[xbin];for(int i =0;i<xbin;i++)rej2[i]=1-eff2[i];

  TGraph *gr2 = new TGraph(xbin,eff1,rej2);
  gr2->SetName("gr2");
  gr2->SetTitle("Z'->ttbar");
  gr2->SetLineColor(kBlue);
  gr2->SetLineWidth(2);
  gr2->GetXaxis()->SetTitle("Eff_{sig}");
  gr2->GetYaxis()->SetTitle("1-Eff_{bkg}");
  gr2->GetXaxis()->SetRangeUser(0,1);
  gr2->GetYaxis()->SetRangeUser(0,1);
  gr2->SetFillStyle(0);

  double rej1[xbin];for(int i =0;i<xbin;i++)rej1[i]=1-eff1[i];

  TGraph *gr3 = new TGraph(xbin,eff1,rej1);
  gr3->SetName("gr3");
  gr3->SetTitle("x=1-y");
  gr3->SetLineColor(kBlue);
  gr3->SetLineStyle(3);
  gr3->SetLineWidth(1);
  gr3->GetXaxis()->SetTitle("Eff_{sig}");
  gr3->GetYaxis()->SetTitle("1-Eff_{bkg}");
  gr3->GetXaxis()->SetRangeUser(0,1);
  gr3->GetYaxis()->SetRangeUser(0,1);
  gr3->SetFillStyle(0);


  TCanvas *cROC1 = new TCanvas("cROC1","cROC1",600,600);

  TMultiGraph *ROC = new TMultiGraph("ROC","ROC");
  ROC->Add(gr2);
  ROC->Add(gr3);

  ROC->Draw("APL");
  ROC->GetXaxis()->SetTitle("Eff_{sig}");
  ROC->GetYaxis()->SetTitle("1-Eff_{bkg}");
  ROC->GetXaxis()->SetLimits(0.,1.);
  ROC->SetMinimum(0.);
  ROC->SetMaximum(1.);
  ROC->SetTitle("ROC, cut: # of btagged Mj>1");

  cROC1->BuildLegend(0.2,0.2,0.2+0.2,0.2+0.15);

  gPad->Update();

  cROC1->SaveAs("plots/cROC1.pdf");


  //chi

  TCanvas *c0 = new TCanvas("c0","c0",600,600);

  h1 = new TH1D("h1","h1",40,-26,-6);
  h2 = new TH1D("h2","h2",40,-26,-6);

  ts->Draw("log(FatJetInfo.Jet_SD_chi)>>h1","FatJetInfo.Jet_SD_chi>0");
  tb->Draw("log(FatJetInfo.Jet_SD_chi)>>h2","FatJetInfo.Jet_SD_chi>0");

  h1->SetLineColor(kBlue);
  h2->SetLineColor(kRed);

  h1->SetTitle("TprimetH M1000");
  h2->SetTitle("Zprimett M1000");

  h1->Scale(1/h1->GetEntries());
  h2->Scale(1/h2->GetEntries());

  h = new THStack("h","");
  h->Add(h1,"sames");
  h->Add(h2,"sames");
  h->Draw("nostack");
  h->GetXaxis()->SetTitle("Log(#chi)");

  leg = new TLegend(0.2,0.65,0.4,0.85);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(h1,"T'->tH M1000","L");
  leg->AddEntry(h2,"Z'->ttbar ","L");
  leg->Draw("SAME");

  gPad->Update();

  tps1 = (TPaveStats*) h1->FindObject("stats");
  tps1->SetTextColor(kBlue);
  X1 = tps1->GetX1NDC();
  Y1 = tps1->GetY1NDC();
  X2 = tps1->GetX2NDC();
  Y2 = tps1->GetY2NDC();
  tps2 = (TPaveStats*) h2->FindObject("stats");
  tps2->SetTextColor(kRed);
  tps2->SetX1NDC(X1);
  tps2->SetX2NDC(X2);
  tps2->SetY1NDC(Y1-(Y2-Y1));
  tps2->SetY2NDC(Y1);
  X1 = tps2->GetX1NDC();
  Y1 = tps2->GetY1NDC();
  X2 = tps2->GetX2NDC();
  Y2 = tps2->GetY2NDC();

  gPad->Update();

  c0->SaveAs("plots/c0.pdf");

  //ROC

  // const int xbin = 40; double xmin = -26; double xmax = -6;

  // double denom1 = 29006; //signal, with chi>0 cut
  // double denom2 = 4003; //bkg, with chi>0 cut

  // double num1[xbin];
  // double num2[xbin];

  // double eff1[xbin];
  // double eff2[xbin];

  // double x[xbin];
  // double inc = (xmax-xmin)/xbin;

  for(int i=0; i<xbin;i++){
    num1[i]=h1->Integral(i+1,xbin+1);
    num2[i]=h2->Integral(i+1,xbin+1);
    eff1[i]=num1[i];
    eff2[i]=num2[i];
    x[i] = xmin + inc*i;
    cout<<i+1<<". x = "<<x[i]<< ", h1 =  "<< num1[i]<< ", h2 = "<< num2[i] << endl;
  }

  for(int i =0;i<xbin;i++)rej2[i]=1-eff2[i];

  gr2 = new TGraph(xbin,eff1,rej2);
  gr2->SetName("gr2");
  gr2->SetTitle("Z'->ttbar");
  gr2->SetLineColor(kBlue);
  gr2->SetLineWidth(2);
  gr2->GetXaxis()->SetTitle("Eff_{sig}");
  gr2->GetYaxis()->SetTitle("1-Eff_{bkg}");
  gr2->GetXaxis()->SetRangeUser(0,1);
  gr2->GetYaxis()->SetRangeUser(0,1);
  gr2->SetFillStyle(0);

  for(int i =0;i<xbin;i++)rej1[i]=1-eff1[i];

  gr3 = new TGraph(xbin,eff1,rej1);
  gr3->SetName("gr3");
  gr3->SetTitle("x=1-y");
  gr3->SetLineColor(kBlue);
  gr3->SetLineStyle(3);
  gr3->SetLineWidth(1);
  gr3->GetXaxis()->SetTitle("Eff_{sig}");
  gr3->GetYaxis()->SetTitle("1-Eff_{bkg}");
  gr3->GetXaxis()->SetRangeUser(0,1);
  gr3->GetYaxis()->SetRangeUser(0,1);
  gr3->SetFillStyle(0);


  TCanvas *cROC = new TCanvas("cROC","cROC",600,600);

  ROC = new TMultiGraph("ROC","ROC");
  ROC->Add(gr2);
  ROC->Add(gr3);

  ROC->Draw("APL");
  ROC->GetXaxis()->SetTitle("Eff_{sig}");
  ROC->GetYaxis()->SetTitle("1-Eff_{bkg}");
  ROC->GetXaxis()->SetLimits(0.,1.);
  ROC->SetMinimum(0.);
  ROC->SetMaximum(1.);
  ROC->SetTitle("ROC Curves");

  cROC->BuildLegend(0.2,0.2,0.2+0.2,0.2+0.15);

  gPad->Update();

  cROC->SaveAs("plots/cROC.pdf");

  //pt

  TCanvas *c2 = new TCanvas("c2","c2",600,600);

  h1 = new TH1D("h1","h1",50,0,1000);
  h2 = new TH1D("h2","h2",50,0,1000);

  ts->Draw("FatJetInfo.Jet_pt>>h1");
  tb->Draw("FatJetInfo.Jet_pt>>h2");

  h1->SetLineColor(kBlue);
  h2->SetLineColor(kRed);

  h1->SetTitle("TprimetH M1000");
  h2->SetTitle("Zprimett M1000");

  h1->Scale(1/h1->GetEntries());
  h2->Scale(1/h2->GetEntries());

  h = new THStack("h","");
  h->Add(h1,"sames");
  h->Add(h2,"sames");
  h->Draw("nostack");
  h->GetXaxis()->SetTitle("FatJet p_{T}");

  leg = new TLegend(0.6,0.65,0.8,0.85);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(h1,"T'->tH M1000","L");
  leg->AddEntry(h2,"Z'->ttbar ","L");
  leg->Draw("SAME");

  gPad->Update();

  tps1 = (TPaveStats*) h1->FindObject("stats");
  tps1->SetTextColor(kBlue);
  X1 = tps1->GetX1NDC();
  Y1 = tps1->GetY1NDC();
  X2 = tps1->GetX2NDC();
  Y2 = tps1->GetY2NDC();
  tps2 = (TPaveStats*) h2->FindObject("stats");
  tps2->SetTextColor(kRed);
  tps2->SetX1NDC(X1);
  tps2->SetX2NDC(X2);
  tps2->SetY1NDC(Y1-(Y2-Y1));
  tps2->SetY2NDC(Y1);
  X1 = tps2->GetX1NDC();
  Y1 = tps2->GetY1NDC();
  X2 = tps2->GetX2NDC();
  Y2 = tps2->GetY2NDC();

  gPad->Update();

  c2->SaveAs("plots/c2.pdf");

  //pt

  TCanvas *c3 = new TCanvas("c3","c3",600,600);

  h1 = new TH1D("h1","h1",50,0,1000);
  h2 = new TH1D("h2","h2",50,0,1000);

  ts->Draw("FatJetInfo.Jet_pt>>h1","FatJetInfo.Jet_SD_chi>0");
  tb->Draw("FatJetInfo.Jet_pt>>h2","FatJetInfo.Jet_SD_chi>0");

  h1->SetLineColor(kBlue);
  h2->SetLineColor(kRed);

  h1->SetTitle("TprimetH M1000");
  h2->SetTitle("Zprimett M1000");

  h1->Scale(1/h1->GetEntries());
  h2->Scale(1/h2->GetEntries());

  h = new THStack("h","");
  h->Add(h1,"sames");
  h->Add(h2,"sames");
  h->Draw("nostack");
  h->GetXaxis()->SetTitle("FatJet p_{T}");
  h->SetTitle("cut : #chi > 0 ");

  leg = new TLegend(0.6,0.65,0.8,0.85);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(h1,"T'->tH M1000","L");
  leg->AddEntry(h2,"Z'->ttbar ","L");
  leg->Draw("SAME");

  gPad->Update();

  tps1 = (TPaveStats*) h1->FindObject("stats");
  tps1->SetTextColor(kBlue);
  X1 = tps1->GetX1NDC();
  Y1 = tps1->GetY1NDC();
  X2 = tps1->GetX2NDC();
  Y2 = tps1->GetY2NDC();
  tps2 = (TPaveStats*) h2->FindObject("stats");
  tps2->SetTextColor(kRed);
  tps2->SetX1NDC(X1);
  tps2->SetX2NDC(X2);
  tps2->SetY1NDC(Y1-(Y2-Y1));
  tps2->SetY2NDC(Y1);
  X1 = tps2->GetX1NDC();
  Y1 = tps2->GetY1NDC();
  X2 = tps2->GetX2NDC();
  Y2 = tps2->GetY2NDC();

  gPad->Update();

  c3->SaveAs("plots/c3.pdf");


}
