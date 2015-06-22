#define Zprime_bkg_eff_cxx
#include "Zprime_bkg_eff.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void Zprime_bkg_eff::Loop()
{
//   In a ROOT session, you can do:
//      root> .L Zprime_bkg_eff.C
//      root> Zprime_bkg_eff t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   const int bin = 50;
   TH1D *h1 = new TH1D("h1","h1",bin,0,1);
   TH1D *h2 = new TH1D("h2","h2",bin,0,1000);
   TH1D *h3 = new TH1D("h3","h3",bin,-23,-7);
   TH1D *h4 = new TH1D("h4","h4",bin,-1,1);

   double denom = 0;
   double num[bin] ={};
   bool pass;
   bool junk;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

   if(jentry%10000==0) cout << "Processed " << jentry
			    << " events" << endl;

      //LOOP over FatJets
      for(int ifj=0; ifj < FatJetInfo_nJet; ifj++){

	junk  = FatJetInfo_Jet_pt[ifj]< 300 ||  FatJetInfo_Jet_massSoftDrop[ifj] > 155 || FatJetInfo_Jet_massSoftDrop[ifj] < 95 ;

	if(junk)continue;
	denom++;
	h1->AddBinContent(0); //underflow bin
	h2->Fill(FatJetInfo_Jet_pt[ifj]);
	h3->Fill(log(FatJetInfo_Jet_SD_chi[ifj]));
	h4->Fill(FatJetInfo_Jet_DoubleSV[ifj]);

	for(int icut=0; icut<bin; icut++){
	    pass = min(PrunedSubJetInfo_Jet_CombIVF[0],PrunedSubJetInfo_Jet_CombIVF[1])>((icut*1.) /(1.*bin));
	    if(pass){
	      num[icut]++;
	      h1->AddBinContent(icut+1);
	    }

	}//end LOOP cut
      }//end LOOP fj

   }//end LOOP event

   cout << "denom = " << denom  << " , h1(0) = " << h1->GetBinContent(0) << endl;
   double x[bin]; for (int i = 0; i<bin;i++){
     x[i]= (1.*i) / (1.*bin);
     cout << "num["<<x[i]<<"] = "<<num[i]<<endl;
   }
   TGraph *gr = new TGraph(bin,x,num);
   //gr->Draw("AP*");

   //h1->Draw();

   TFile* f1 = new TFile("Zprime_bkg_eff.root","recreate");
   gr->Write();

   h1->GetXaxis()->SetTitle("double CombIVF eff");
   h2->GetXaxis()->SetTitle("Fatjet p_{T}");
   h3->GetXaxis()->SetTitle("SD Log( #chi)");
   h4->GetXaxis()->SetTitle("DoubleSV eff");

   //   h2->Draw();

   h1->Write();
   h2->Write();
   h3->Write();
   h4->Write();

   f1->Write();
   f1->Close();


}
