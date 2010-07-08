{
	TFile f1("hist_lnQ_b1235_sys.root");
	TFile f2("hist_lnQ_nosys.root");
	TH1D *h1 = (TH1D*)f1->Get("h");
	TH1D *h2 = (TH1D*)f2->Get("h");
	TCanvas c;
	c.SetLogy(1);
	h1->SetStats(0);
	h2->SetStats(0);
	h1->Rebin(1000);
	h2->Rebin(1000);
	h1->SetNormFactor(1);
	h2->SetNormFactor(1);
	h1->GetXaxis()->SetRangeUser(-60, 60);
	h2->GetXaxis()->SetRangeUser(-60, 60);
	h2->Draw();
	h2->SetLineColor(kRed);
	h1->Draw("same");
	h2->SetTitle(";-2lnQ (b-only hypothesis); Entries (normalized to 1)");

	TFile f3("likelihood_pdf_sys.root");
	TFile f4("likelihood_pdf.root");
	TCanvas *c3 = (TCanvas*)f3->Get("c");
	TCanvas *c4 = (TCanvas*)f4->Get("c");
	TGraph *g3=(TGraph*)c3->GetListOfPrimitives()->FindObject("Graph");
	TGraph *g4=(TGraph*)c4->GetListOfPrimitives()->FindObject("Graph");
	TCanvas c5;
	g3->Draw("AL");
	g4->SetLineColor(kRed);
	g4->Draw("Lsame");

}

