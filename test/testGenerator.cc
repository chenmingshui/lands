#include <iostream>
#include "TPad.h"
#include "TPaveStats.h"
#include "CLsLimit.h"
#include "CRandom.h"
#include "PlotUtilities.h"
//#include "FloridaStyle.C"
#include <TError.h>
#include "TSystem.h"
#include "CountingModel.h"
#include "BayesianBase.h"
#include "LimitBands.h"
using std::cout;
using std::endl;
using namespace lands;
int main(int argc, const char* argv[]){
	double s=1.;
	double b=0.46;
	double d=0.2;
	double s_err=0;
	double b_err=1.;
	int EsEb_correlated=1;

	CRandom *rdm = new CRandom(1234);  //initilize a random generator
	CountingModel *cms=new CountingModel();
	cms->SetRdm(rdm);
	cms->AddChannel(s,b); 
	cms->AddObservedData(0,d); 
	cms->AddUncertainty(0,              0,            s_err,         1,                       1                ); 
	cms->AddUncertainty(0,              1,            b_err,         1,                       EsEb_correlated==0?2:1        ); 
	if(s_err !=0 || b_err !=0 ) { 
		cms->SetUseSystematicErrors(true);
	}

	CountingModel *cmsG=new CountingModel();
	cmsG->SetRdm(rdm);
	cmsG->AddChannel(s,b); 
	cmsG->AddObservedData(0,d); 
	cmsG->AddUncertainty(0,              0,            s_err,         2,                       1                ); 
	cmsG->AddUncertainty(0,              1,            b_err,         2,                       EsEb_correlated==0?2:1        ); 
	if(s_err !=0 || b_err !=0 ) { 
		cmsG->SetUseSystematicErrors(true);
	}


	BayesianBase bys(cms, 0.05, 1.e-3);
	bys.SetNumToys(20000);
	bys.SetDebug(0);
	double rtmp;
	rtmp = bys.Limit();
	cout<<"------------------------------------------------------------"<<endl;
	cout<<"Observed Upper Limit on the ratio R at 95\% CL = "<<rtmp<<endl;
	cout<<"------------------------------------------------------------"<<endl;

	// get the background distribution by throwing 1M toys
	TH1F *hb=new TH1F("LogNormal",";b; entries",1000, 0, b+7*b*b_err);	
	TH1F *hbG=new TH1F("TrunGaus",";b; entries",1000, 0, b+7*b*b_err);	
	TH1F *hboverb=new TH1F("logNormal",";#upsilon=b/b_{0}; entries",1000, 0, 5);	
	TH1F *hbGoverb=new TH1F("trunGaus",";#upsilon=b/b_{0}; entries",1000, 0, 5);	
	TH1F *hlnb=new TH1F("lnb",";ln(b); entries",1000, -10, 10);	
	vector< vector<double> > vv;
	for(int i=0; i<1000000; i++){
		vv = cms->FluctuatedNumbers();
		hb->Fill(vv[0][1]);
		hboverb->Fill(vv[0][1]/b);
		hlnb->Fill(log(vv[0][1]));
		vv = cmsG->FluctuatedNumbers();
		hbG->Fill(vv[0][1]);
		hbGoverb->Fill(vv[0][1]/b);
	}
	TCanvas c;
	c.SetLogx(0);
	hb->SetLineColor(kRed);
	hbG->SetLineColor(kBlue);
	hb->Draw();
	hbG->Draw("sames");
	gPad->Update();
	TPaveStats *st=(TPaveStats*)gPad->GetPrimitive("stats");
	st->SetName("hbGstat");
	st->SetY1NDC(0.650);
	st->SetY2NDC(0.80);
//	st->SetTextColor(kBlue);
	c.Modified();


	Save(&c,"b");
	c.SetLogx(1);
	Save(&c,"b_logx");

	c.SetLogx(0);
	hlnb->Draw();
	Save(&c,"lnb");


	hboverb->SetLineColor(kRed);
	hbGoverb->SetLineColor(kBlue);
	hboverb->Draw();
	hbGoverb->Draw("sames");
	gPad->Update();
	st=(TPaveStats*)gPad->GetPrimitive("stats");
	st->SetName("hbGstat");
	st->SetY1NDC(0.650);
	st->SetY2NDC(0.80);
//	st->SetTextColor(kBlue);
	c.SetLogy(0);
	c.SetLogx(0);
	Save(&c,"rate");
	c.SetLogx(1);
	Save(&c,"rate_logx");
	c.SetLogx(0);
	c.SetLogy(1);
	Save(&c,"rate_logy");

	bys.PosteriorPdf();
	TPaveText *pt = SetTPaveText(0.2, 0.7, 0.3, 0.9);
	DrawEvolution2D pdfr(bys.GetVR(), bys.GetVP(), ";r;Likelihood", "likelihood_pdf", pt);	
	pdfr.draw();


	return 1;
}
