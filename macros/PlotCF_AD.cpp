TH1F* SetMyStyle(TH1F *const &hist){
  hist->GetXaxis()->SetTitleSize(0.055);
  hist->GetYaxis()->SetTitleSize(0.055);
  hist->GetXaxis()->SetTitleOffset(0.7);
  hist->GetYaxis()->SetTitleOffset(0.75);
  hist->GetXaxis()->SetNdivisions(208);
  hist->GetYaxis()->SetNdivisions(208);
  hist->SetMarkerStyle(kFullCircle);
  hist->SetMarkerSize(1.);
  hist->SetMarkerColor(kRed+2);
  hist->SetLineColor(kRed+2);
  return hist;
}

void PlotCF_AD() {
  TFile *fi = new TFile("Pi.root","READ");
  Int_t markerColor[] = {kBlack,kRed+2,kGreen+2,kBlue+2};
  Int_t markerStyle[] = {kFullCircle, kOpenSquare, kOpenTriangleUp, kOpenTriangleDown};
  const Int_t nCentBin = 4;
  const TString CentralityBinning[nCentBin] = {"No AD cut","AD>4cm","AD>6cm","AD>8cm"};
  TH1F *tmp[nCentBin], *tmp2[nCentBin];
  for (Int_t icent(0);icent<nCentBin;icent++) {
    if (icent == 0) {
      tmp[icent] = (TH1F*) fi->Get(Form("hQinvNum_SL_0_1"));
      tmp2[icent] = (TH1F*) fi->Get(Form("hQinvDen_SL_0_1"));
      // tmp[icent] = (TH1F*) fi->Get(Form("hQinvNum_0"));
      // tmp2[icent] = (TH1F*) fi->Get(Form("hQinvDen_0"));
    } else {
      tmp[icent] = (TH1F*) fi->Get(Form("hQinvNum_SL_AD_0_%i",icent-1));
      tmp2[icent] = (TH1F*) fi->Get(Form("hQinvDen_SL_AD_0_%i",icent-1));
      // tmp[icent] = (TH1F*) fi->Get(Form("hQinvNum_AD_0_%i",icent-1));
      // tmp2[icent] = (TH1F*) fi->Get(Form("hQinvDen_AD_0_%i",icent-1));
    }
    tmp[icent]->Scale(1./tmp[icent]->Integral());
    tmp2[icent]->Scale(1./tmp2[icent]->Integral());
    tmp[icent]->Divide(tmp2[icent]);
    tmp[icent]->GetXaxis()->SetTitle("q_{inv}, GeV/c");
    tmp[icent]->GetYaxis()->SetTitle("C(q_{inv})");
    tmp[icent]->GetXaxis()->SetRangeUser(0,0.2);
    tmp[icent]->GetYaxis()->SetRangeUser(0.,1.4);
    SetMyStyle(tmp[icent]);
    tmp[icent]->SetMarkerColor(markerColor[icent]);
    tmp[icent]->SetLineColor(markerColor[icent]);
    tmp[icent]->SetMarkerStyle(markerStyle[icent]);
    
    // cout << "123" << endl;
  }
  gStyle->SetOptStat(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetErrorX(0);
  TCanvas c;
  tmp[0]->Draw("L");
  tmp[1]->Draw("L same");
  tmp[2]->Draw("L same");
  tmp[3]->Draw("L same");
  TLine one;
  one.SetLineStyle(2);
  one.DrawLine(0,1,0.5,1);
  TLatex tex;
  tex.DrawLatexNDC(0.3,0.8,Form("STAR Au+Au #sqrt{s_{NN}} = 27 GeV, 0-10%%"));
  tex.SetTextFont(42);
  tex.DrawLatexNDC(0.3,0.7,Form("#pi^{+}#pi^{+}, 0.2 < k_{T} < 1.0 GeV/c"));
  TLegend *leg = new TLegend(0.3,0.2,0.6,0.4);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.04);
  for (Int_t icent(0);icent<nCentBin;icent++) leg->AddEntry(tmp[icent],CentralityBinning[icent].Data(), "p");
  leg->Draw();
  c.SaveAs(Form("CF_PiPi_AD_cent_0_10.pdf"));
  
}