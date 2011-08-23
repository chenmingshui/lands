#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TArrow.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TString.h>
#include <TError.h>
#include <TStyle.h>
#include <math.h>
#include <iostream>
#include "TMath.h"
#include <algorithm>

#include "PlotUtilities_standalone.h"
using namespace std;
//------------------------------------------------------------------------------
// SetTPaveText
//------------------------------------------------------------------------------
TPaveText* SetTPaveText(Double_t x1, Double_t y1, Double_t x2, Double_t y2)
{
	TPaveText* t = new TPaveText(x1, y1, x2, y2, "ndc");

	t->SetBorderSize(    0);
	t->SetFillStyle (    0);
	t->SetTextAlign (   12);
	t->SetTextFont  (   42);
	t->SetTextSize  (0.035);
	//  t->SetTextColor (kRed);

	return t;
}
//------------------------------------------------------------------------------
// DrawLegend
//------------------------------------------------------------------------------
TLegend* DrawLegend(Float_t x1,
		Float_t y1,
		TH1F*   hist,
		TString label,
		TString option,
		Float_t tsize,
		Float_t xoffset,
		Float_t yoffset)
{
	TLegend* legend = new TLegend(x1,
			y1,
			x1 + xoffset,
			y1 + yoffset);

	legend->SetBorderSize(    0);
	legend->SetFillColor (    0);
	legend->SetFillStyle (    0);//transparent
	legend->SetTextAlign (   12);
	legend->SetTextFont  (   42);
	legend->SetTextSize  (tsize);

	legend->AddEntry(hist, label.Data(), option.Data());
	legend->Draw();

	return legend;
}
TLegend* DrawLegend(Float_t x1,
		Float_t y1,
		TGraph*   hist,
		TString label,
		TString option,
		Float_t tsize,
		Float_t xoffset,
		Float_t yoffset)
{
	TLegend* legend = new TLegend(x1,
			y1,
			x1 + xoffset,
			y1 + yoffset);

	legend->SetBorderSize(    0);
	legend->SetFillColor (    0);
	legend->SetFillStyle (    0);//transparent
	legend->SetTextAlign (   12);
	legend->SetTextFont  (   42);
	legend->SetTextSize  (tsize);

	legend->AddEntry(hist, label.Data(), option.Data());
	legend->Draw();

	return legend;
}
void MoveLegend(TLegend *leg, Float_t x1, 
	Float_t y1, Float_t xoffset, Float_t yoffset){
	leg->SetX1(x1);
	leg->SetX2(x1+xoffset);
	leg->SetY1(y1);
	leg->SetY2(y1+yoffset);
}
void DrawText(double x1, double y1, double x2, double y2, TString text, double size ){
	TPaveText *pt = SetTPaveText(x1, y1, x2, y2);
	pt->AddText(text);
	pt->SetTextSize(size);
	pt->Draw();
}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++              reading file with format
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
/*
   TGraph::TGraph(const char *filename, const char *format, Option_t *)
   : TNamed("Graph",filename), TAttLine(), TAttFill(1,1001), TAttMarker()
   {
// Graph constructor reading input from filename
// filename is assumed to contain at least two columns of numbers
// the string format is by default "%lg %lg"

Double_t x,y;
ifstream infile(filename);
if(!infile.good()){
MakeZombie();
Error("TGraph", "Cannot open file: %s, TGraph is Zombie",filename);
fNpoints = 0;
} else {
fNpoints = 100;  //initial number of points
}
if (!CtorAllocate()) return;
std::string line;
Int_t np=0;
while(std::getline(infile,line,'\n')){
if(2 != sscanf(line.c_str(),format,&x,&y) ) {
continue; // skip empty and ill-formed lines
}
SetPoint(np,x,y);
np++;
}
Set(np);
}

 */
PlotWithBelts::PlotWithBelts(
		double *r1sigLt, double *r1sigHt, 
		double *r2sigLt, double *r2sigHt, 
		double *rmeant, int npointst,  
		double *xpointst,  string  ssavet, TPaveText *ptt,
		double xmint,  double xmaxt, double ymint, double ymaxt, bool logYt,
		string stitlet){
	//double robs[1000];

	xmin = xmint;
	xmax = xmaxt;
	ymin = ymint;
	ymax = ymaxt;
	ssave = ssavet;
	logY = logYt;
	pt=ptt;

	npoints=npointst;
	r1sigL=r1sigLt;
	r1sigH=r1sigHt;
	r2sigL=r2sigLt;
	r2sigH=r2sigHt;
	rmean=rmeant;
	xpoints=xpointst;
	stitle=stitlet;

	bDrawGrobs=2; // don't draw
}
PlotWithBelts::PlotWithBelts(
		double *r1sigLt, double *r1sigHt, 
		double *r2sigLt, double *r2sigHt, 
		double *rmeant, double *robst, int npointst,  
		double *xpointst,  string  ssavetmp, TPaveText *pttmp,
		double xmintmp,  double xmaxtmp, double ymintmp, double ymaxtmp, bool logYtmp,
		string stitlet){

	xmin = xmintmp;
	xmax = xmaxtmp;
	ymin = ymintmp;
	ymax = ymaxtmp;
	ssave = ssavetmp;
	logY = logYtmp;
	pt=pttmp;

	npoints=npointst;
	r1sigL=r1sigLt;
	r1sigH=r1sigHt;
	r2sigL=r2sigLt;
	r2sigH=r2sigHt;
	rmean=rmeant;
	robs=robst;
	xpoints=xpointst;
	stitle=stitlet;


	bDrawGrobs=1;

}
void PlotWithBelts::drawLegend(string s_gr, string s_gGreen, string s_gYellow, string s_grobs){
	delete legend;
	legend = DrawLegend(0.2, 0.3, gr, s_gr.c_str(), "l");
	legend->AddEntry(gGreen,s_gGreen.c_str(), "f");
	legend->AddEntry(gYellow,s_gYellow.c_str(), "f");
	if(bDrawGrobs==1){
		legend->AddEntry(grobs, s_grobs.c_str(), "lp");
	}
}
PlotWithBelts::~PlotWithBelts(){
	clear();
}
void PlotWithBelts::save(){
	string seps = ssave+".eps";
	string sgif = ssave+".gif";
	string sroot = ssave+".root";
	cCanvas->Print(sroot.c_str());
	cCanvas->Print(seps.c_str());
	cCanvas->Print(sgif.c_str());
}
void PlotWithBelts::clear(){

	if(gr)	delete gr; gr=NULL;
	if(gYellow)	delete gYellow;
	delete gGreen;
	//	delete pt;
	//	if(lineOne)	delete lineOne;
	if(grobs)	delete grobs;
	if(legend)	delete legend;

	//	if(ptCMSPreli)	delete ptCMSPreli;
	//	pt=NULL;
	//	cCanvas=NULL;
	/*
	   gr=NULL;
	   gYellow=NULL;
	   grobs=NULL;
	   legend=NULL;
	   lineOne=NULL;
	 */
	lineOne=NULL;
	if(hframe) delete hframe;
}
void PlotWithBelts::plot(){
	cCanvas = new TCanvas(ssave.c_str(),"Canvas");
	//cCanvas->DrawFrame(xmin,ymin,xmax, ymax, stitle.c_str());
	TString stmp = "hframe"; stmp += ssave;
	hframe= new TH1F(stmp, stitle.c_str(), 1000, xmin, xmax);
	hframe->SetMinimum(ymin);
	hframe->SetMaximum(ymax);
	hframe->SetStats(0);
	hframe->SetFillStyle(1);
	hframe->Draw(" ");
	
	cCanvas->SetLogy(logY);


	gr = new TGraph(npoints, xpoints, rmean); 
	gYellow = new TGraph(2*npoints);
	for(int n=0; n<npoints; n++){
		gYellow->SetPoint(n, xpoints[n], r2sigH[n]);
		gYellow->SetPoint(npoints+n, xpoints[npoints-n-1], r2sigL[npoints-n-1]);
	}
	gYellow->SetFillColor(kYellow);
	gYellow->SetLineColor(kYellow);
	gYellow->Draw("f");

	gGreen = new TGraph(2*npoints);
	for(int n=0; n<npoints; n++){
		gGreen->SetPoint(n, xpoints[n], r1sigH[n]);
		gGreen->SetPoint(npoints+n, xpoints[npoints-n-1], r1sigL[npoints-n-1]);
	}
	gGreen->SetFillColor(kGreen);
	gGreen->SetLineColor(kGreen);
	gGreen->Draw("f");

	lineOne = new TLine(xmin,1, xmax, 1);
	lineOne->SetLineWidth(2);
	lineOne->SetLineStyle(1);
	lineOne->SetLineColor(kBlack);
	lineOne->Draw("same");

	gr->SetLineWidth(1);
	gr->SetLineStyle(2);
	gr->Draw("LP");


	grobs=NULL;
	//	if(bDrawGrobs!=2) bDrawGrobs=1; // if it has been assigned in another constructor, then, keep it as it is. 
	if(bDrawGrobs==1){ // need NOTE re-code here
		grobs=new TGraph(npoints, xpoints, robs);
		grobs->SetMarkerStyle(21);
		grobs->SetMarkerColor(kBlue);
		grobs->SetLineWidth(2);
		grobs->SetLineColor(kBlue);
		grobs->Draw("LP");
	}

	legend = DrawLegend(0.2, 0.3, gr, "bkgd-only: mean", "l");
	legend->AddEntry(gGreen,"bkgd-only: 1 #sigma band", "f");
	legend->AddEntry(gYellow,"bkgd-only: 2 #sigma band", "f");
	if(bDrawGrobs==1){
		legend->AddEntry(grobs,"observed", "lp");
	}

	DrawCMS();
	pt->Draw();
	/*ptCMSPreli=SetTPaveText(0.7, 0.9, 0.8, 1.0);
	  ptCMSPreli->AddText("CMS Preliminary");
	  ptCMSPreli->SetTextSize(0.05);
	  ptCMSPreli->Draw();*/

//	save();

}
DrawSigBkgPdfs::DrawSigBkgPdfs(double (*pdfs)(double *, double *), double *pars, int npars,
		double (*pdfb)(double *, double *), double *parb, int nparb, 
		double Ns, double Nb,
		double xstart, double xstop, bool logY, string ssave, string stitle, bool debug){
	_pdfs=pdfs;	_pdfb=pdfb; 	_pars=pars;	_parb=parb;	_npars=npars;	_nparb=nparb;
	_xstart=xstart;	_xstop=xstop;	_Ns=Ns; _Nb=Nb; _logY=logY; _ssave=ssave; _stitle=stitle; _debug=debug;
	fs=0;fb=0;hs=0;hb=0;htot=0;cCanvas=0;legend=0; 
}
DrawSigBkgPdfs::~DrawSigBkgPdfs(){
	delete fs; delete fb; delete hs; delete hb; delete htot; delete legend;
	cCanvas=NULL; 
}
void DrawSigBkgPdfs::draw(){
	cCanvas=new TCanvas("c","c");
	cCanvas->cd();
	cCanvas->SetLogy(_logY);
	fs=new TF1("fs",_pdfs, _xstart, _xstop, _npars);
	fs->SetParameters(_pars);
	fb=new TF1("fb",_pdfb, _xstart, _xstop, _nparb);
	fb->SetParameters(_parb);
	fs->SetNpx(1000);
	fb->SetNpx(1000);
	if(_debug){
		fs->Print();
		cout<<"fs integral = "<<fs->Integral(_xstart,_xstop)<<endl;
		fb->Print();
		cout<<"fb integral = "<<fb->Integral(_xstart,_xstop)<<endl;
	}
	hs = (TH1F*)fs->GetHistogram();
	hb = (TH1F*)fb->GetHistogram();
	hb->Scale(_Nb/(_Ns+_Nb));
	htot=(TH1F*)hs->Clone("htot");//fs->GetHistogram();
	htot->Scale(_Ns/(_Ns+_Nb));
	htot->Add(hb);
	htot->SetTitle(_stitle.c_str());
	htot->Draw();
	hs->Scale(_Ns/(_Ns+_Nb));
	if(_debug)cout<<"hist s+b integral = "<<htot->Integral("width")<<endl;
	if(_debug)cout<<"hist s   integral = "<<hs->Integral("width")<<endl;
	if(_debug)cout<<"hist b   integral = "<<hb->Integral("width")<<endl;
	hs->SetLineColor(kBlue);
	hb->SetLineColor(kRed);
	hs->Draw("same");
	hb->Draw("same");

	legend = DrawLegend(0.2, 0.3, htot, "sig+bkg", "l");
	legend->AddEntry(hs,"sig", "l");
	legend->AddEntry(hb,"bkg", "l");

//	save();
}
void DrawSigBkgPdfs::save(){
	string seps = _ssave+".eps";
	string sgif = _ssave+".gif";
	string sroot = _ssave+".root";
	cCanvas->Print(sroot.c_str());
	cCanvas->Print(seps.c_str());
	cCanvas->Print(sgif.c_str());
}
void DrawSigBkgPdfs::drawLegend(string s_tot, string s_s, string s_b){
	delete legend;
	legend = DrawLegend(0.2, 0.3, htot, s_tot.c_str(), "l");
	legend->AddEntry(hs,s_s.c_str(), "l");
	legend->AddEntry(hb,s_b.c_str(), "l");
}

DrawEvolution2D::DrawEvolution2D(vector<double> vx, vector<double> vy, string stitle, string ssave, TPaveText *pt, bool debug){
	_vx=vx; _vy=vy; _stitle=stitle; _ssave=ssave; _pt=pt; _debug=debug;	
	_logY=0; cCanvas=0; legend=0; lineOne=0; graph=0;
}
DrawEvolution2D::~DrawEvolution2D(){
	delete legend; delete lineOne; delete graph;
	cCanvas=0; _pt=0;
}
void DrawEvolution2D::draw(){
	cCanvas= new TCanvas("c","c");	
	cCanvas->SetLogy(_logY);
	int nr = _vx.size();
	double *rtmp=new double[nr];
	double *cls_btmp = new double[nr];
	for(int i=0; i<nr; i++){
		rtmp[i]=_vx[i]; cls_btmp[i]=_vy[i];
	}
	int *ir=new int[nr];
	TMath::Sort(nr, rtmp, ir, 0);
	double *r=new double[nr]; 
	double *cls_b = new double[nr];
	for(int i=0; i<nr; i++){
		r[i]=rtmp[ir[i]];
		cls_b[i]=cls_btmp[ir[i]];
	}

	graph = new TGraph(nr, r, cls_b);
	graph->SetMarkerStyle(21);
	graph->SetMarkerColor(kBlue);
	graph->SetLineWidth(2);
	graph->SetLineColor(kBlue);
	graph->SetTitle(_stitle.c_str());
	graph->Draw("ALP");

	double xmin=graph->GetXaxis()->GetXmin();
	double xmax=graph->GetXaxis()->GetXmax();
	lineOne = new TLine(xmin, 0.05, xmax, 0.05);
	lineOne->SetLineWidth(2);
	lineOne->SetLineStyle(1);
	lineOne->SetLineColor(kBlack);
	lineOne->Draw("same");

	_pt->Draw();

//	save();

	delete [] rtmp; 
	delete [] cls_btmp;
	delete [] ir;
	delete [] r; 
	delete [] cls_b;
}
void DrawEvolution2D::save(){
	Save(cCanvas,_ssave);
	/*string seps = _ssave+".eps";
	  string sgif = _ssave+".gif";
	  string sroot = _ssave+".root";
	  cCanvas->Print(sroot.c_str());
	  cCanvas->Print(seps.c_str());
	  cCanvas->Print(sgif.c_str());*/
}
PlotXvsCummulativeProb::PlotXvsCummulativeProb(
		double* vx_not_sorted, int nexps,
		string ssave, string stitle, TPaveText *pt){
	vx_not_sorted =0 ;
	//SortAndCumulative(vx_not_sorted, nexps, _vx, _vy);
	//GetBandsByLinearInterpolation(_vx, _vy, _1SigmaLow,_1SigmaHigh, _2SigmaLow, _2SigmaHigh);
	_ssave=ssave; _stitle=stitle;_pt=pt;
	legend=0;lineOne=0;graph=0;cCanvas=0;diffHist=0;bShowErrorBar=true; _nexps=nexps;
	_graphOption="";
}
PlotXvsCummulativeProb::PlotXvsCummulativeProb(
		vector<double> vx_not_sorted,
		string ssave, string stitle, TPaveText *pt){
	vx_not_sorted.clear();
	//SortAndCumulative(vx_not_sorted, _vx, _vy);
	//GetBandsByLinearInterpolation(_vx, _vy, _1SigmaLow,_1SigmaHigh, _2SigmaLow, _2SigmaHigh);
	_ssave=ssave; _stitle=stitle;_pt=pt;
	legend=0;lineOne=0;graph=0;cCanvas=0;diffHist=0; bShowErrorBar=true; _nexps=vx_not_sorted.size();
	_graphOption="";
}
PlotXvsCummulativeProb::PlotXvsCummulativeProb(
		vector<double> vx_sorted,
		vector<double> vcumulaP_sorted,
		double m1s, double p1s, double m2s, double p2s,
		string ssave, string stitle, TPaveText *pt){
	_vx=vx_sorted; _vy=vcumulaP_sorted; _1SigmaLow=m1s; _1SigmaHigh=p1s; _2SigmaLow=m2s; _2SigmaHigh=p2s;
	_ssave=ssave; _stitle=stitle;_pt=pt;
	legend=0;lineOne=0;graph=0;cCanvas=0;diffHist=0;bShowErrorBar=false;
	_graphOption="";
}
PlotXvsCummulativeProb::~PlotXvsCummulativeProb(){
	delete legend; delete lineOne; delete graph; cCanvas=0; _pt=0; delete diffHist;
}
void PlotXvsCummulativeProb::draw(){
	// ---- need to check everything is ok... 
	// ---- or re-sort it no matter it's done 
	cCanvas=new TCanvas("c1","c1");
	int nsize=(int)_vy.size();
	double xmax= _vx[nsize-1]+1;
	cCanvas->DrawFrame(0,0,xmax,1, _stitle.c_str());

	TBox *boxGreen  = new TBox(_1SigmaLow,0, _1SigmaHigh,1);
	TBox *boxYellow = new TBox(_2SigmaLow,0, _2SigmaHigh,1);
	boxYellow->SetFillColor(kYellow);
	boxYellow->Draw("same");
	boxGreen->SetFillColor(kGreen);
	boxGreen->Draw("same");

	double GreenBandLow = (1- 0.683)/2.; //1 sigma
	double GreenBandHigh = 1 - (1- 0.683)/2.; //1 sigma
	double YellowBandLow = (1- 0.955)/2.; //2 sigma
	double YellowBandHigh = 1 - (1- 0.955)/2.; //2 sigma
	//double median = 0.5;

	TLine *lineGrnLow = new TLine(0,GreenBandLow, xmax, GreenBandLow);
	TLine *lineGrnHigh = new TLine(0,GreenBandHigh, xmax, GreenBandHigh);
	TLine *lineYlwLow = new TLine(0,YellowBandLow, xmax, YellowBandLow);
	TLine *lineYlwHigh = new TLine(0,YellowBandHigh, xmax, YellowBandHigh);
	TLine *lineMedian= new TLine(0,0.5, xmax, 0.5);
	lineGrnLow->SetLineColor(kGreen+1);
	lineGrnLow->SetLineWidth(2);
	lineGrnLow->Draw("same");
	lineGrnHigh->SetLineColor(kGreen+1);
	lineGrnHigh->SetLineWidth(2);
	lineGrnHigh->Draw("same");
	lineYlwLow->SetLineColor(kYellow+1);
	lineYlwLow->SetLineWidth(2);
	lineYlwLow->Draw("same");
	lineYlwHigh->SetLineColor(kYellow+1);
	lineYlwHigh->SetLineWidth(2);
	lineYlwHigh->Draw("same");
	lineMedian->SetLineColor(kBlack);
	lineMedian->SetLineWidth(2);
	lineMedian->Draw("same");

	double *rn=new double[nsize];
	double *pn=new double[nsize];
	double *er=new double[nsize];
	double *ep=new double[nsize];


	//	double rn[1000000], pn[1000000];	
	//	double er[1000000], ep[1000000]; //cause segment fault
	for(int i=0; i<nsize; i++){
		rn[i]=_vx[i];
		pn[i]=_vy[i];
		er[i]=0;ep[i]=0;
		if(bShowErrorBar)ep[i]=sqrt(pn[i]*(1-pn[i])/(double)_nexps);
	}
	graph = new TGraphErrors(nsize,rn,pn,er, ep);
	graph->SetMarkerStyle(8);
	graph->SetMarkerColor(kRed);
	graph->SetMarkerSize(0.3);
	if(_graphOption != ""){ graph->Draw(_graphOption.c_str());
	}
	else{
		if(!bShowErrorBar)graph->Draw("CP");
		else graph->Draw("LP"); //if no A, then no title displayed
	}

	_pt->Draw();
	DrawCMS();
//	save();

	delete [] rn;
	delete [] pn;
	delete [] er;
	delete [] ep;

}
void PlotXvsCummulativeProb::save(){
	Save(cCanvas,_ssave);
}
void PlotXvsCummulativeProb::drawDifferencial(string stitle, string ssave){
	int nsize=(int)_vy.size();

	double rn[10000], pn[10000];	
	double xmax= _vx[nsize-1]+1;
	diffHist=new TH1F("h",stitle.c_str(), 100, 0, xmax+1);
	diffHist->SetStats(0);
	for(int i=0; i<nsize; i++){
		rn[i]=_vx[i];
		if(i>0){
			pn[i]=_vy[i]*_nexps - _vy[i-1]*_nexps;
		}else{
			pn[i]=_vy[i]*_nexps;
		}
		for(int j=0; j<(int)pn[i]; j++)
			diffHist->Fill(rn[i]);
	}
	cCanvas=new TCanvas("c22","c22");
	diffHist->Draw("E");
	TPaveText *stat=SetTPaveText(0.7, 0.8, 0.8, 0.9);
	char cmean[256];
	sprintf(cmean, "Mean = %.2f #pm %.2f", diffHist->GetMean(), diffHist->GetMeanError());
	char crms[256];
	sprintf(crms,  "RMS  = %.2f ", diffHist->GetRMS());
	stat->AddText(cmean);
	stat->AddText(crms);
	stat->Draw();


	_pt->Draw();

	DrawCMS();
	Save(cCanvas, ssave);
}

void DrawCMS(string cms, double x1, double y1, double x2, double y2, double txtsize){
	cms = "";
	x1 = x2;
	y1 = y2;
	txtsize=1;
	// cms="CMS Preliminary";
	/*
	   TPaveText *ptCMSPreli=SetTPaveText(x1, y1, x2, y2);
	   ptCMSPreli->AddText(cms.c_str());
	   ptCMSPreli->SetTextSize(txtsize);
	   ptCMSPreli->Draw();
	 */
}
void Save(TCanvas *cCanvas, string _ssave){
	string seps = _ssave+".eps";
	string spdf = _ssave+".pdf";
	string sgif = _ssave+".gif";
	string sroot = _ssave+".root";
	cCanvas->Print(sroot.c_str());
	cCanvas->Print(seps.c_str());
	cCanvas->Print(spdf.c_str());
	cCanvas->Print(sgif.c_str());
}
DrawPdfM2logQ::DrawPdfM2logQ(vector<double> vPdfM2logQ_sb, vector<double> vPdfM2logQ_b, 
		double m2lnQ_d, string sdata, string stitle, string ssave, TPaveText *pt){ // show the famous -2lnQ with s, b, d
	_vsb=vPdfM2logQ_sb; _vb=vPdfM2logQ_b; _m2lnQ_d=m2lnQ_d; _lineTitle=sdata; _stitle=stitle; _ssave=ssave; _pt=pt;
	cCanvas=0;legend=0;lineOne=0;gSB=0;gB=0;hSB=0;hB=0;
	bFromPseudoExps=true;
}
DrawPdfM2logQ::DrawPdfM2logQ(vector< pair<double,double> > vPdfM2logQ_sb, vector< pair<double,double> > vPdfM2logQ_b, 
		double m2lnQ_d, string sdata, string stitle, string ssave, TPaveText *pt){ // show the famous -2lnQ with s, b, d
	_vpsb=vPdfM2logQ_sb; _vpb=vPdfM2logQ_b; _m2lnQ_d=m2lnQ_d; _lineTitle=sdata; _stitle=stitle; _ssave=ssave; _pt=pt;
	cCanvas=0;legend=0;lineOne=0;gSB=0;gB=0;hSB=0;hB=0;
	bFromPseudoExps=false;
}
DrawPdfM2logQ::~DrawPdfM2logQ(){
	cCanvas=0; delete legend; delete lineOne; delete gSB; delete gB; delete hSB; delete hB; _pt=0;
}
void DrawPdfM2logQ::draw(){
	if(bFromPseudoExps) draw2Hist();
	else draw2Graph();
}
void DrawPdfM2logQ::draw2Hist(){
	cCanvas = new TCanvas("cM2logQ","cM2logQ");

	TH1F *hPdfM2logQ= new TH1F("hPdfM2logQ","",100, 0, 0);
	for(int i=0; i<(int)_vsb.size(); i++){
		hPdfM2logQ->Fill(_vsb[i]);
		hPdfM2logQ->Fill(_vb[i]);
	}
	hSB = new TH1F("hPdfM2logQ_sb","",100, hPdfM2logQ->GetXaxis()->GetXmin(), hPdfM2logQ->GetXaxis()->GetXmax());
	hB  = new TH1F("hPdfM2logQ_b", "",100, hPdfM2logQ->GetXaxis()->GetXmin(), hPdfM2logQ->GetXaxis()->GetXmax());
	delete hPdfM2logQ;
	hSB->SetStats(0);
	hB->SetStats(0);
	for(int i=0; i<(int)_vsb.size(); i++){
		hSB->Fill(_vsb[i]);
		hB->Fill(_vb[i]);
	}

	hSB->SetLineColor(kBlue);
	hSB->SetLineStyle(3);
	hB->SetLineColor(kRed);
	hB->SetTitle(_stitle.c_str());
	hB->Draw();
	hSB->Draw("same");

	double ymin, ymax;
	ymin=1e-2;
	ymax=hSB->GetMaximum();
	if(ymax<hB->GetMaximum()) ymax=hB->GetMaximum();
	lineOne = new TLine(_m2lnQ_d, ymin, _m2lnQ_d, ymax);
	lineOne->SetLineWidth(2);
	lineOne->SetLineStyle(1);
	lineOne->SetLineColor(kBlack);
	lineOne->Draw("same");

	legend = DrawLegend(0.2, 0.3, hSB, "s+b", "lf");
	legend->AddEntry(hB, "b-only", "lf");
	legend->AddEntry(lineOne, _lineTitle.c_str(), "l");

	_pt->Draw();
	Save(cCanvas, _ssave);
}
void DrawPdfM2logQ::draw2Graph(){
	int nq_sb = _vpsb.size();
	int nq_b = _vpb.size();

	double *q_sb=new double[nq_sb];
	double *q_b=new double[nq_b];
	double *p_sb=new double[nq_sb];
	double *p_b=new double[nq_b];
	for(int i=0; i<nq_sb; i++) {q_sb[i]=_vpsb[i].first; p_sb[i]=_vpsb[i].second;}
	for(int i=0; i<nq_b; i++) {q_b[i]=_vpb[i].first; p_b[i]=_vpb[i].second;}

	cCanvas= new TCanvas("cM2logQ","cM2logQ");
	gSB=new TGraph(nq_sb, q_sb, p_sb);
	gB=new TGraph(nq_b, q_b, p_b);

	// --- first draw a frame 
	double xmin, xmax, ymin,ymax;
	xmin=gSB->GetXaxis()->GetXmin();
	if(xmin>gB->GetXaxis()->GetXmin()) xmin=gB->GetXaxis()->GetXmin();
	xmax=gSB->GetXaxis()->GetXmax();
	if(xmax<gB->GetXaxis()->GetXmax()) xmax=gB->GetXaxis()->GetXmax();
	ymin=gSB->GetYaxis()->GetXmin();
	if(ymin>gB->GetYaxis()->GetXmin()) ymin=gB->GetYaxis()->GetXmin();
	ymax=gSB->GetYaxis()->GetXmax();
	if(ymax<gB->GetYaxis()->GetXmax()) ymax=gB->GetYaxis()->GetXmax();
	double x[4]={xmin,xmin,xmax,xmax};
	double y[4]={ymin,ymax,ymin,ymax};
	TGraph *frame=new TGraph(4,x,y);	
	frame->SetTitle(_stitle.c_str());
	frame->Draw("AP");
	//cout<<" xmin="<<xmin<<" xmax="<<xmax<<" ymin="<<ymin<<" ymax="<<ymax<<endl;

	gSB->SetLineColor(kBlue);
	gSB->SetLineStyle(3);
	gSB->SetMarkerStyle(24);
	gSB->SetMarkerColor(kBlue);
	gB->SetLineColor(kRed);
	gB->SetMarkerStyle(20);
	gB->SetMarkerColor(kRed);
	gB->SetTitle(_stitle.c_str());
	gB->Draw("CP");
	gSB->Draw("CP");

	lineOne = new TLine(_m2lnQ_d, ymin, _m2lnQ_d, ymax);
	lineOne->SetLineWidth(2);
	lineOne->SetLineStyle(1);
	lineOne->SetLineColor(kBlack);
	lineOne->Draw("same");

	legend = DrawLegend(0.2, 0.25, gSB, "s+b", "lp");
	legend->AddEntry(gB, "b-only", "lp");
	legend->AddEntry(lineOne, _lineTitle.c_str(), "l");

	_pt->Draw();
	Save(cCanvas, _ssave);

	delete frame; // if delete this frame, then we may not be able to modify it later
	delete [] q_sb;
	delete [] q_b;
	delete [] p_sb;
	delete [] p_b;
}
void DrawPdfM2logQ::drawLegend(string s_sb, string s_b, string s_d){
	delete legend;
	if(bFromPseudoExps){
		legend = DrawLegend(0.2, 0.3, hSB, s_sb.c_str(), "lf");
		legend->AddEntry(hB,s_b.c_str(), "lf");
		legend->AddEntry(lineOne,s_d.c_str(), "l");
	}
	else{
		legend = DrawLegend(0.2, 0.25, gSB, s_sb.c_str(), "lp");
		legend->AddEntry(gB, s_b.c_str(), "lp");
		legend->AddEntry(lineOne, s_d.c_str(), "l");
	}
}
DrawMultiGraph::DrawMultiGraph(double *x1, double *y1, int n1, string title1, 
		string stitle, string ssave, TPaveText *pt,
		double xmin, double xmax, double ymin, double ymax, bool logY){
	_x[0]=x1;_y[0]=y1;_n[0]=n1;_title[0]=title1;
	_stitle=stitle;_ssave=ssave;_pt=pt;
	_xmin=xmin; _xmax=xmax; _ymin=ymin; _ymax=ymax; _logY=logY;
	_numGraphAdded=1; cCanvas=0; legend=NULL;lineOne=0;
	for(int i=0; i<maxNumGraphs; i++) { 
		_gr[i]=0;
		_drawOptions[i]="";
	}
	legendX=0.2; legendY=0.3;
}
DrawMultiGraph::DrawMultiGraph(double *x1, double *y1, int n1, string title1, double *x2, double *y2, int n2, string title2, 
		string stitle, string ssave, TPaveText *pt,
		double xmin, double xmax, double ymin, double ymax, bool logY){
	_x[0]=x1;_y[0]=y1;_n[0]=n1;_title[0]=title1;
	_x[1]=x2;_y[1]=y2;_n[1]=n2;_title[1]=title2;
	_stitle=stitle;_ssave=ssave;_pt=pt;
	_xmin=xmin; _xmax=xmax; _ymin=ymin; _ymax=ymax; _logY=logY;
	_numGraphAdded=2; cCanvas=0; legend=NULL;lineOne=0;
	for(int i=0; i<maxNumGraphs; i++) { 
		_gr[i]=0;
		_drawOptions[i]="";
	}
	legendX=0.2; legendY=0.3;
}
DrawMultiGraph::~DrawMultiGraph(){
	for(int i=0; i<maxNumGraphs; i++){
		if(_gr[i])	delete _gr[i];
	}
	if(lineOne)	lineOne=0;
	if(legend)  delete legend; 
	if(cCanvas) cCanvas=0; 
	if(_pt)_pt=0;
}
void DrawMultiGraph::add(double *x, double *y, int n, string title){
	if(_numGraphAdded>=maxNumGraphs) {cout<<"Error: too many graphs added, do nothing"<<endl; return;}
	_numGraphAdded++;
	_x[_numGraphAdded-1]=x; _y[_numGraphAdded-1]=y;  _n[_numGraphAdded-1]=n; _title[_numGraphAdded-1]=title;	
}
void DrawMultiGraph::draw(){
	cCanvas=new TCanvas("cExcLimits","cExcLimits");
	cCanvas->DrawFrame(_xmin,_ymin,_xmax, _ymax, _stitle.c_str());
	cCanvas->SetLogy(_logY);
	for(int i=0; i<_numGraphAdded; i++){
		_gr[i]=new TGraph(_n[i], _x[i], _y[i] );
		_gr[i]->SetLineWidth(1);
		_gr[i]->SetLineStyle(_glineStyle[i]);
		_gr[i]->SetLineColor(_glineColor[i]);
		_gr[i]->SetMarkerStyle(_gmarkerStyle[i]);
		_gr[i]->SetMarkerColor(_glineColor[i]);
		if(_drawOptions[i]!="")_gr[i]->Draw(_drawOptions[i].c_str());
		else _gr[i]->Draw("LP");
	}


	if(_drawOptions[0]!=""){
		legend = DrawLegend(legendX, legendY, _gr[0], _title[0].c_str(), _drawOptions[0].c_str());
	}else{
		legend = DrawLegend(legendX, legendY, _gr[0], _title[0].c_str(), "lp");
	}

	for(int i=1; i<_numGraphAdded; i++) {
		if(_drawOptions[i]!=""){
			legend->AddEntry(_gr[i], _title[i].c_str(), _drawOptions[i].c_str());
		}else{
			legend->AddEntry(_gr[i], _title[i].c_str(), "lp");
		}
	}

	//	cout<<"delete me legendY2="<<legend->GetY2()<<endl;

	lineOne = new TLine(_xmin,1, _xmax, 1);
	lineOne->SetLineWidth(2);
	lineOne->SetLineStyle(1);
	lineOne->SetLineColor(kBlack);
	lineOne->Draw("same");

	_pt->Draw();

	DrawCMS();
	//	Save(cCanvas,_ssave); // it changes the legend 
	//	cout<<"delete me legendY2="<<legend->GetY2()<<endl;
}
DrawPdfRLikelihood::DrawPdfRLikelihood(double r95, double *par, int npar, string stitle, string ssave, TPaveText *pt){
	_r95=r95; _par=par; _npar=npar; _stitle=stitle; _ssave=ssave; _pt=pt;
	hist=0; cCanvas=0; arrow=0;  _bUnbinned=false;
}
DrawPdfRLikelihood::~DrawPdfRLikelihood(){
	delete hist; cCanvas=0; delete arrow; _pt=0;
}
void DrawPdfRLikelihood::draw(){
	double xmax = _r95*2;
	double step = xmax/100.;
	double y[100];
	double ymax=0;
	double norm=0;
	hist=new TH1F("h",_stitle.c_str(),100,0,xmax);
	for(int i=0; i<100; i++){
		//double x[1]={step*i};		
		//if(!_bUnbinned)y[i]=P_n0_Given_brs(x, _par, _npar);
		if(ymax<y[i]) ymax=y[i];
		norm+=y[i]*step;
		hist->SetBinContent(i,y[i]);
	}	
	ymax/=norm;
	hist->Scale(1./norm);
	hist->SetStats(0);// don't show the lands box
	cCanvas=new TCanvas("c1","c1");
	hist->Draw("C");	

	arrow = new TArrow(_r95, ymax/7.,_r95,0,0.04);
	arrow->SetLineWidth(3);
	arrow->SetLineColor(kRed);
	arrow->Draw();

	_pt->Draw();

	DrawCMS();
	Save(cCanvas,_ssave);
}

void GetLimits(TTree *tree, vector<double>& inputMH, vector<double>& inputLimits, vector<double> & inputLimitErrs){	
   // Declaration of leaf types
   Double_t        mH;
   Double_t        limit;
   Double_t        limitErr;
   Double_t        significance;
   Double_t        pvalue;
   Double_t        rm2s;
   Double_t        rm1s;
   Double_t        rmedian;
   Double_t        rmean;
   Double_t        rp1s;
   Double_t        rp2s;

   // List of branches
   TBranch        *b_mH;   //!
   TBranch        *b_limit;   //!
   TBranch        *b_limitErr;   //!
   TBranch        *b_significance;   //!
   TBranch        *b_pvalue;   //!
   TBranch        *b_rm2s;   //!
   TBranch        *b_rm1s;   //!
   TBranch        *b_rmedian;   //!
   TBranch        *b_rmean;   //!
   TBranch        *b_rp1s;   //!
   TBranch        *b_rp2s;   //!

   TTree *fChain = tree;
   if(tree->GetBranch("mH")){
	   fChain->SetBranchAddress("mH", &mH, &b_mH);
   }
   if(tree->GetBranch("mh")){
	   fChain->SetBranchAddress("mh", &mH, &b_mH);
   }

   fChain->SetBranchAddress("limit", &limit, &b_limit);
   fChain->SetBranchAddress("limitErr", &limitErr, &b_limitErr);
   fChain->SetBranchAddress("significance", &significance, &b_significance);
   fChain->SetBranchAddress("pvalue", &pvalue, &b_pvalue);
   fChain->SetBranchAddress("rm2s", &rm2s, &b_rm2s);
   fChain->SetBranchAddress("rm1s", &rm1s, &b_rm1s);
   fChain->SetBranchAddress("rmedian", &rmedian, &b_rmedian);
   fChain->SetBranchAddress("rmean", &rmean, &b_rmean);
   fChain->SetBranchAddress("rp1s", &rp1s, &b_rp1s);
   fChain->SetBranchAddress("rp2s", &rp2s, &b_rp2s);

   Long64_t nentries = tree->GetEntries();


   vector< vector<double> > vv_sameMH; vv_sameMH.clear();
   //Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
	   Long64_t ientry = tree->GetEntry(jentry);
	   if (ientry < 0) break;
	   vector<double>::iterator id = std::find(inputMH.begin(), inputMH.end(), mH); 
	   if(id!=inputMH.end()){
		   vv_sameMH[id-inputMH.begin()].push_back(limit);
	   }else{
		   vector<double> v; v.clear();
		   v.push_back(limit);vv_sameMH.push_back(v);
		   inputLimits.push_back(limit);
		   inputLimitErrs.push_back(limitErr);
		   inputMH.push_back(mH);
	   }
   }
   for(int i=0; i<vv_sameMH.size();i++){ 
	   if(vv_sameMH[i].size()>1) {
		   std::sort(vv_sameMH[i].begin(), vv_sameMH[i].end());
		   double avr=0, avrerr=0;
		   for(int j=0; j<vv_sameMH[i].size(); j++){
			   avr+=vv_sameMH[i][j];
		   }
		   avr/=(float)(vv_sameMH[i].size());

		   for(int j=0; j<vv_sameMH[i].size(); j++){
			   avrerr+=(vv_sameMH[i][j]-avr)*(vv_sameMH[i][j]-avr);
		   }
		   avrerr = TMath::Sqrt(avrerr)/TMath::Sqrt( vv_sameMH[i].size() * (vv_sameMH[i].size()-1) );
		   //inputLimits[i]=avr;
		   inputLimits[i]=vv_sameMH[i][(int)(vv_sameMH[i].size()*0.5)];
		   inputLimitErrs[i]=avrerr;
	   }
   }
}
void GetPValues(TTree *tree, vector<double>& inputMH, vector<double>& inputLimits, vector<double> & inputLimitErrs_m2s, vector<double> & inputLimitErrs_m1s, vector<double>&inputLimitErrs_p1s, vector<double>&inputLimitErrs_p2s ){
	// Declaration of leaf types
	Double_t        mH;
	Double_t        limit;
	Double_t        limitErr;
	Double_t        significance;
	Double_t        pvalue;

	// List of branches
	TBranch        *b_mH;   //!
	TBranch        *b_limit;   //!
	TBranch        *b_limitErr;   //!
	TBranch        *b_significance;   //!
	TBranch        *b_pvalue;   //!

	TTree *fChain = tree;
	if(tree->GetBranch("mH")){
		fChain->SetBranchAddress("mH", &mH, &b_mH);
		fChain->SetBranchAddress("pvalue", &pvalue, &b_pvalue);
	}
	if(tree->GetBranch("mh")){
		fChain->SetBranchAddress("mh", &mH, &b_mH);
		fChain->SetBranchAddress("limit", &pvalue, &b_limit);
	}

	//   fChain->SetBranchAddress("limit", &limit, &b_limit);
	//   fChain->SetBranchAddress("limitErr", &limitErr, &b_limitErr);
	//   fChain->SetBranchAddress("significance", &significance, &b_significance);

	Long64_t nentries = tree->GetEntries();


   vector< vector<double> > vv_sameMH; vv_sameMH.clear();
   //Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
	   Long64_t ientry = tree->GetEntry(jentry);
	   if (ientry < 0) break;
	   vector<double>::iterator id = std::find(inputMH.begin(), inputMH.end(), mH); 
	   if(id!=inputMH.end()){
		   vv_sameMH[id-inputMH.begin()].push_back(pvalue);
	   }else{
		   vector<double> v; v.clear();
		   v.push_back(pvalue);vv_sameMH.push_back(v);
		   inputLimits.push_back(pvalue);
		   inputLimitErrs_m2s.push_back(pvalue);
		   inputLimitErrs_m1s.push_back(pvalue);
		   inputLimitErrs_p2s.push_back(pvalue);
		   inputLimitErrs_p1s.push_back(pvalue);
		   inputMH.push_back(mH);
	   }
   }
   for(int i=0; i<vv_sameMH.size();i++){ 
	   if(vv_sameMH[i].size()>1) {
		   cout<<" More than 1 pvalues at this mass "<<inputMH[i]<<endl;
		   exit(1);
		   std::sort(vv_sameMH[i].begin(), vv_sameMH[i].end());
		   double avr=0, avrerr=0;
		   for(int j=0; j<vv_sameMH[i].size(); j++){
			   avr+=vv_sameMH[i][j];
		   }
		   avr/=(float)(vv_sameMH[i].size());

		   for(int j=0; j<vv_sameMH[i].size(); j++){
			   avrerr+=(vv_sameMH[i][j]-avr)*(vv_sameMH[i][j]-avr);
		   }
		   avrerr = TMath::Sqrt(avrerr)/TMath::Sqrt( vv_sameMH[i].size() * (vv_sameMH[i].size()-1) );
		   //inputLimits[i]=avr;
		   inputLimits[i]=vv_sameMH[i][(int)(vv_sameMH[i].size()*0.5)];
		   inputLimitErrs_m2s[i]=vv_sameMH[i][(int)(vv_sameMH[i].size()*0.0275)];
		   inputLimitErrs_m1s[i]=vv_sameMH[i][(int)(vv_sameMH[i].size()*0.16)];
		   inputLimitErrs_p2s[i]=vv_sameMH[i][(int)(vv_sameMH[i].size()*0.975)];
		   inputLimitErrs_p1s[i]=vv_sameMH[i][(int)(vv_sameMH[i].size()*0.84)];
		   //inputLimitErrs[i]=avrerr;
	   }
   }
}

void GetLimitBands(TTree *tree, vector<double>& inputMH, vector<double>& inputLimits, vector<double> & inputLimitsM2S, vector<double>& inputLimitsM1S, 
		vector<double>& inputLimitsMEDIAN, vector<double>& inputLimitsP1S, vector<double>& inputLimitsP2S){	
   // Declaration of leaf types
   Double_t        mH;
   Double_t        limit;
   Double_t        limitErr;
   Double_t        significance;
   Double_t        pvalue;
   Double_t        rm2s;
   Double_t        rm1s;
   Double_t        rmedian;
   Double_t        rmean;
   Double_t        rp1s;
   Double_t        rp2s;

   // List of branches
   TBranch        *b_mH;   //!
   TBranch        *b_limit;   //!
   TBranch        *b_limitErr;   //!
   TBranch        *b_significance;   //!
   TBranch        *b_pvalue;   //!
   TBranch        *b_rm2s;   //!
   TBranch        *b_rm1s;   //!
   TBranch        *b_rmedian;   //!
   TBranch        *b_rmean;   //!
   TBranch        *b_rp1s;   //!
   TBranch        *b_rp2s;   //!

   TTree *fChain = tree;
   if(tree->GetBranch("mH")){
	   fChain->SetBranchAddress("mH", &mH, &b_mH);
   }
   if(tree->GetBranch("mh")){
	   fChain->SetBranchAddress("mh", &mH, &b_mH);
   }

   fChain->SetBranchAddress("limit", &limit, &b_limit);
   fChain->SetBranchAddress("limitErr", &limitErr, &b_limitErr);
   fChain->SetBranchAddress("significance", &significance, &b_significance);
   fChain->SetBranchAddress("pvalue", &pvalue, &b_pvalue);
   fChain->SetBranchAddress("rm2s", &rm2s, &b_rm2s);
   fChain->SetBranchAddress("rm1s", &rm1s, &b_rm1s);
   fChain->SetBranchAddress("rmedian", &rmedian, &b_rmedian);
   fChain->SetBranchAddress("rmean", &rmean, &b_rmean);
   fChain->SetBranchAddress("rp1s", &rp1s, &b_rp1s);
   fChain->SetBranchAddress("rp2s", &rp2s, &b_rp2s);

   Long64_t nentries = tree->GetEntries();

   //Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
	   Long64_t ientry = tree->GetEntry(jentry);
	   if (ientry < 0) break;
	   inputLimits.push_back(limit);
	   inputMH.push_back(mH);
	   inputLimitsM2S.push_back(rm2s);
	   inputLimitsM1S.push_back(rm1s);
	   inputLimitsP2S.push_back(rp2s);
	   inputLimitsP1S.push_back(rp1s);
	   inputLimitsMEDIAN.push_back(rmedian);
   }
}
void GetMuHat(TTree *tree, vector<double>& inputMH, vector<double>& inputLimits ){
	// Declaration of leaf types
	Double_t        mH;
	Double_t        limit;

	// List of branches
	TBranch        *b_mH;   //!
	TBranch        *b_limit;   //!

	TTree *fChain = tree;
	if(tree->GetBranch("mH")){
		fChain->SetBranchAddress("mH", &mH, &b_mH);
		fChain->SetBranchAddress("rmean", &limit, &b_limit);
	}
	if(tree->GetBranch("mh")){
		fChain->SetBranchAddress("mh", &mH, &b_mH);
		fChain->SetBranchAddress("limit", &limit, &b_limit);
	}


	Long64_t nentries = tree->GetEntries();

	for (Long64_t jentry=0; jentry<nentries;jentry++) {
		Long64_t ientry = tree->GetEntry(jentry);
		if (ientry < 0) break;
		inputLimits.push_back(limit);
		inputMH.push_back(mH);
	}
}
