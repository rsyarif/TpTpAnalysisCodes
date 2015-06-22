{
  TFile *fs = new TFile("Tprime_sig_eff_matched.root");
  //TFile *fb = new TFile("Zprime_bkg_eff.root");
  TFile *fb = new TFile("TT_Mtt_bkg_eff_unmatched.root");

  const int  bin = 50;

  gStyle->SetOptStat("");

 //      _             _     _         _____                _    _______      ________
 //     | |           | |   | |       / ____|              | |  |_   _\ \    / /  ____|
 //   __| | ___  _   _| |__ | | ___  | |     ___  _ __ ___ | |__  | |  \ \  / /| |__
 //  / _` |/ _ \| | | | '_ \| |/ _ \ | |    / _ \| '_ ` _ \| '_ \ | |   \ \/ / |  __|
 // | (_| | (_) | |_| | |_) | |  __/ | |___| (_) | | | | | | |_) || |_   \  /  | |
 //  \__,_|\___/ \__,_|_.__/|_|\___|  \_____\___/|_| |_| |_|_.__/_____|   \/   |_|


  TH1D *hs = (TH1D*)fs->Get("h1");
  TH1D *hb = (TH1D*)fb->Get("h1");

  //const int bin = hs->GetSize() - 2 ;

  double x[bin] ;//= {};
  double y[bin] ;//= {};

  double denom_s = hs->GetBinContent(0) * 1.;
  double denom_b = hb->GetBinContent(0) * 1.;

  for(int i = 0; i<bin; i++){
    x[i] = (hs->GetBinContent(i+1) * 1.) / denom_s ; //cout << "x["<<i<<"] = "<<x[i]<<endl;
    y[i] = 1 - ( (hb->GetBinContent(i+1) * 1.) / denom_b ) ; //cout << "y["<<i<<"] = "<<y[i]<<endl;
  }

  TGraph *gr1 = new TGraph(bin,x,y);
  gr1->SetName("gr1");
  gr1->SetTitle("double CombIVF");
  gr1->SetLineColor(kBlue);
  gr1->SetLineWidth(2);
  gr1->GetXaxis()->SetTitle("Eff_{sig}");
  gr1->GetYaxis()->SetTitle("1-Eff_{bkg}");
  gr1->GetXaxis()->SetLimits(0,1);
  gr1->SetMinimum(0.);
  gr1->SetMaximum(1.);
  gr1->SetFillStyle(0);
  //gr1->Draw("ALP");
  //c1->SaveAs("plots/ROC_doubleSubjetCSVIVF.pdf");

  TCanvas *cs1 = new TCanvas("cs1","cs1",600,600);
  cs1->cd();

  hs->SetLineColor(kBlue);
  hb->SetLineColor(kRed);

  hs->Scale(1/hs->GetBinContent(0));
  hb->Scale(1/hb->GetBinContent(0));

  THStack *h = new THStack("h","");
  h->Add(hs,"sames");
  h->Add(hb,"sames");
  h->Draw("nostack");
  h->GetXaxis()->SetTitle("double CombIVF");
  h->SetTitle("p_{T}>300, 95<m_{SoftDrop}<155");

  TLegend *leg = new TLegend(0.2,0.65,0.4,0.85);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(hs,"T'->tH M1000","L");
  //leg->AddEntry(hb,"Z'ttbar M1000","L");
  leg->AddEntry(hb,"TT_Mtt M1000","L");
  leg->Draw("SAME");

  gPad->Update();

  cs1->SaveAs("plots/doubleCombIVF_sigbkg.eps");

 //            _          __   _____  _______      ____        _____  _______      _____   __
 //    	      (_)        / /  / ____|/ ____\ \    / /_ |      / ____|/ ____\ \    / /__ \  \ \
 //  _ __ _ _  _ _ __   | |  | |    | (___  \ \  / / | |     | |    | (___  \ \  / /   ) |  | |
 // | '_ ` _ \| | '_ \  | |  | |     \___ \  \ \/ /  | |     | |     \___ \  \ \/ /   / /   | |
 // | | | | | | | | | | | |  | |____ ____) |  \  /   | |  _  | |____ ____) |  \  /   / /_   | |
 // |_| |_| |_|_|_| |_| | |   \_____|_____/    \/    |_| ( )  \_____|_____/    \/   |____|  | |
 //                      \_\                             |/                                /_/

  TH1D *hs1_1 = (TH1D*)fs->Get("h1_1");
  TH1D *hb1_1 = (TH1D*)fb->Get("h1_1");

  //const int bin = hs1_1->GetSize() - 2 ;

  double x1_1[bin] ;//= {};
  double y1_1[bin] ;//= {};

  double denom_s1_1 = hs1_1->GetEntries() * 1.;
  double denom_b1_1 = hb1_1->GetEntries() * 1.;

  for(int i = 0; i<bin; i++){
    x1_1[i] = (hs1_1->GetBinContent(i+1) * 1.) / denom_s1_1 ; //cout << "x1_1["<<i<<"] = "<<x1_1[i]<<endl;
    y1_1[i] = 1 - ( (hb1_1->GetBinContent(i+1) * 1.) / denom_b1_1 ) ; //cout << "y1_1["<<i<<"] = "<<y1_1[i]<<endl;
  }

  TGraph *gr1_1 = new TGraph(bin,x1_1,y1_1);
  gr1_1->SetName("gr1_1");
  gr1_1->SetTitle("min(CSV_{sj1},CSV_{sj2})");
  gr1_1->SetLineColor(kBlue);
  gr1_1->SetLineWidth(2);
  gr1_1->GetXaxis()->SetTitle("Eff_{sig}");
  gr1_1->GetYaxis()->SetTitle("1-Eff_{bkg}");
  gr1_1->GetXaxis()->SetLimits(0,1);
  gr1_1->SetMinimum(0.);
  gr1_1->SetMaximum(1.);
  gr1_1->SetFillStyle(0);
  //gr1_1->Draw("ALP");
  //c1->SaveAs("plots/ROC_doubleSubjetCSVIVF.pdf");

  TCanvas *cs1_1 = new TCanvas("cs1_1","cs1_1",600,600);
  cs1_1->cd();

  hs1_1->SetLineColor(kBlue);
  hb1_1->SetLineColor(kRed);

  hs1_1->Scale(1/hs1_1->GetBinContent(0));
  hb1_1->Scale(1/hb1_1->GetBinContent(0));

  THStack *h1_1 = new THStack("h1_1","");
  h1_1->Add(hs1_1,"sames");
  h1_1->Add(hb1_1,"sames");
  h1_1->Draw("nostack");
  h1_1->GetXaxis()->SetTitle("min(CSV_sj1,CSV_sj2)");
  h1_1->SetTitle("p_{T}>300, 95<m_{SoftDrop}<155");

  TLegend *leg1_1 = new TLegend(0.2,0.65,0.4,0.85);
  leg1_1->SetFillStyle(0);
  leg1_1->SetBorderSize(0);
  leg1_1->AddEntry(hs1_1,"T'->tH M1000","L");
  //leg1_1->AddEntry(hb1_1,"Z'ttbar M1000","L");
  leg1_1->AddEntry(hb1_1,"TT_Mtt M1000","L");
  leg1_1->Draw("SAME");

  gPad->Update();

  cs1_1->SaveAs("plots/minCSV1CSV2_sigbkg.eps");

 //   _____ _____
 //  / ____|  __ \
 // | (___ | |  | |
 //  \___ \| |  | |
 //  ____) | |__| |
 // |_____/|_____/


  TH1D *hs3 = (TH1D*)fs->Get("h3");
  TH1D *hb3 = (TH1D*)fb->Get("h3");

  //const int bin = hs->GetSize() - 2 ;

  //const int  bin = 50;

  double x3[bin] ;//= {};
  double y3[bin] ;//= {};

  double denom_s3 = hs3->GetEntries() * 1.;
  double denom_b3 = hb3->GetEntries() * 1.;

  for(int i = 0; i<bin; i++){
    x3[i] = (hs3->Integral(i+1,bin) * 1.) / denom_s3 ; //cout << "x3["<<i<<"] = "<<x3[i]<< ", " << (hs3->Integral(i+1,bin) * 1.) << " / "<<  denom_s3 << endl;
    y3[i] = 1 - ( (hb3->Integral(i+1,bin) * 1.) / denom_b3 ) ; //cout << "y3["<<i<<"] = "<<y3[i]<< ", " << (hb3->Integral(i+1,bin) * 1.) <<" / " << denom_b3 <<endl;
  }

  TGraph *gr3 = new TGraph(bin,x3,y3);
  gr3->SetName("gr3");
  gr3->SetTitle("SD Hbb tagger");
  gr3->SetLineColor(kGreen+2);
  gr3->SetLineWidth(2);
  gr3->GetXaxis()->SetTitle("Eff_{sig}");
  gr3->GetYaxis()->SetTitle("1-Eff_{bkg}");
  gr3->GetXaxis()->SetLimits(0,1);
  gr3->SetMinimum(0.);
  gr3->SetMaximum(1.);
  gr3->SetFillStyle(0);
  //gr3->Draw("ALP");

  TCanvas *cs3 = new TCanvas("cs3","cs3",600,600);
  cs3->cd();

  hs3->SetLineColor(kBlue);
  hb3->SetLineColor(kRed);

  hs3->Scale(1/hs3->GetEntries());
  hb3->Scale(1/hb3->GetEntries());

  THStack *h3 = new THStack("h3","");
  h3->Add(hs3,"sames");
  h3->Add(hb3,"sames");
  h3->Draw("nostack");
  h3->GetXaxis()->SetTitle("SD log(#chi)");
  h3->SetTitle("p_{T}>300, 95<m_{SoftDrop}<155");

  TLegend *leg3 = new TLegend(0.2,0.65,0.4,0.85);
  leg3->SetFillStyle(0);
  leg3->SetBorderSize(0);
  leg3->AddEntry(hs,"T'->tH M1000","L");
  //leg3->AddEntry(hb,"Z'ttbar M1000","L");
  leg3->AddEntry(hb,"TT_Mtt M1000","L");
  leg3->Draw("SAME");

  gPad->Update();

  cs3->SaveAs("plots/SD_sigbkg.eps");



 //                        _    _ _     _       _
 //                       | |  | | |   | |     | |
 //  _ __   _____      __ | |__| | |__ | |__   | |_ __ _  __ _  __ _  ___ _ __
 // | '_ \ / _ \ \ /\ / / |  __  | '_ \| '_ \  | __/ _` |/ _` |/ _` |/ _ \ '__|
 // | | | |  __/\ V  V /  | |  | | |_) | |_) | | || (_| | (_| | (_| |  __/ |
 // |_| |_|\___| \_/\_/   |_|  |_|_.__/|_.__/   \__\__,_|\__, |\__, |\___|_|
 //                                                       __/ | __/ |
 //                                                      |___/ |___/


  TH1D *hs4 = (TH1D*)fs->Get("h4");
  TH1D *hb4 = (TH1D*)fb->Get("h4");

  //const int bin = hs->GetSize() - 2 ;

  //const int  bin = 50;

  double x4[bin] ;//= {};
  double y4[bin] ;//= {};

  double denom_s4 = hs4->GetEntries() * 1.;
  double denom_b4 = hb4->GetEntries() * 1.;

  for(int i = 0; i<bin; i++){
    x4[i] = (denom_s4 - hs4->Integral(i+1,bin) * 1.) / denom_s4 ; //cout << "x4["<<i<<"] = "<<x4[i]<< ", " << (hs4->Integral(i+1,bin) * 1.) << " / "<<  denom_s4 << endl;
    y4[i] = 1 - ( (denom_b4 - hb4->Integral(i+1,bin) * 1.) / denom_b4 ) ; //cout << "y4["<<i<<"] = "<<y4[i]<< ", " << (hb4->Integral(i+1,bin) * 1.) <<" / " << denom_b4 <<endl;
  }

  TGraph *gr4 = new TGraph(bin,x4,y4);
  gr4->SetName("gr4");
  gr4->SetTitle("new Hbb tagger");
  gr4->SetLineColor(kRed);
  gr4->SetLineWidth(2);
  gr4->GetXaxis()->SetTitle("Eff_{sig}");
  gr4->GetYaxis()->SetTitle("1-Eff_{bkg}");
  gr4->GetXaxis()->SetLimits(0,1);
  gr4->SetMinimum(0.);
  gr4->SetMaximum(1.);
  gr4->SetFillStyle(0);
  //gr4->Draw("ALP");

  TCanvas *cs4 = new TCanvas("cs4","cs4",600,600);
  cs4->cd();

  hs4->SetLineColor(kBlue);
  hb4->SetLineColor(kRed);

  hs4->Scale(1/hs4->GetEntries());
  hb4->Scale(1/hb4->GetEntries());

  THStack *h4 = new THStack("h4","");
  h4->Add(hs4,"sames");
  h4->Add(hb4,"sames");
  h4->Draw("nostack");
  h4->GetXaxis()->SetTitle("new HBB Tagger");
  h4->SetTitle("p_{T}>300, 95<m_{SoftDrop}<155");

  TLegend *leg4 = new TLegend(0.2,0.65,0.4,0.85);
  leg4->SetFillStyle(0);
  leg4->SetBorderSize(0);
  leg4->AddEntry(hs,"T'->tH M1000","L");
  //leg4->AddEntry(hb,"Z'ttbar M1000","L");
  leg4->AddEntry(hb,"TT_Mtt M1000","L");
  leg4->Draw("SAME");

  gPad->Update();

  cs4->SaveAs("plots/newHbbTagger_sigbkg.eps");


 //           _      _
 //     /\   | |    | |
 //    /  \  | |    | |
 //   / /\ \ | |    | |
 //  / ____ \| |____| |____
 // /_/    \_\______|______|

  TCanvas *c1 = new TCanvas("c1","c1",600,600);
  c1->cd();
  //  c1->SetLogy();

  TMultiGraph *ROC = new TMultiGraph("ROC","ROC");
  ROC->Add(gr1);
  ROC->Add(gr3);
  ROC->Add(gr4);
  ROC->Draw("APL");
  ROC->GetXaxis()->SetTitle("Eff_{sig}");
  ROC->GetYaxis()->SetTitle("1-Eff_{bkg}");
  ROC->GetXaxis()->SetLimits(0.,1.);
  ROC->SetMinimum(0.);
  ROC->SetMaximum(1.);
  ROC->SetTitle("ROC");

  c1->BuildLegend(0.2,0.2,0.2+0.2,0.2+0.15);

  gPad->Update();

  c1->SaveAs("plots/ROC_all.eps");


 //     _______
 //    |__   __|
 //  _ __ | |
 // | '_ \| |
 // | |_) | |
 // | .__/|_|
 // | |
 // |_|


  TH1D *hs2 = (TH1D*)fs->Get("h2");
  TH1D *hb2 = (TH1D*)fb->Get("h2");

  TCanvas *cs2 = new TCanvas("cs2","cs2",600,600);
  cs2->cd();

  hs2->SetLineColor(kBlue);
  hb2->SetLineColor(kRed);

  hs2->Scale(1/hs2->GetEntries());
  hb2->Scale(1/hb2->GetEntries());

  THStack *h2 = new THStack("h2","");
  h2->Add(hs2,"sames");
  h2->Add(hb2,"sames");
  h2->Draw("nostack");
  h2->GetXaxis()->SetTitle("Fatjet p_{T}");
  h2->SetTitle("p_{T}>300, 95<m_{SoftDrop}<155");

  TLegend *leg2 = new TLegend(0.2,0.65,0.4,0.85);
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  leg2->AddEntry(hs2,"T'->tH M1000","L");
  //leg2->AddEntry(hb2,"Z'ttbar M1000","L");
  leg2->AddEntry(hb2,"TT_Mtt M1000","L");
  leg2->Draw("SAME");

  gPad->Update();

  cs2->SaveAs("plots/pT_sigbkg.eps");


}
