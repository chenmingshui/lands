#ifndef  PLOTUTILITIES_H
#define  PLOTUTILITIES_H
#include <TString.h>
#include <string>
#include <TPaveText.h>
#include <TLegend.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TArrow.h>

using namespace std;

namespace lands{
void DrawCMS(string cms="CMS Preliminary", double x1=0.7, double y1=0.9, double x2=0.8, double y2=1.0, double txtsize=0.05); //FIXME: don't do this inside class
TPaveText* SetTPaveText       (Double_t           x1,
		Double_t           y1,
		Double_t           x2,
		Double_t           y2);
TLegend*   DrawLegend         (Float_t            x1,
		Float_t            y1,
		TH1F*              hist,
		TString            label,
		TString            option,
		Float_t            tsize   = 0.04,
		Float_t            xoffset = 0.25,
		Float_t            yoffset = 0.15);
TLegend*   DrawLegend         (Float_t            x1,
		Float_t            y1,
		TGraph*              hist,
		TString            label,
		TString            option,
		Float_t            tsize   = 0.04,
		Float_t            xoffset = 0.25,
		Float_t            yoffset = 0.15);
void MoveLegend(TLegend *leg, Float_t x1, 
	Float_t y1, Float_t xoffset=0.25, Float_t yoffset=0.15);
void Save(TCanvas *cCanvas, string _ssave);
void DrawText(double x1, double y1, double x2, double y2, TString text, double size = 0.045);
class PlotWithBelts{
	public:
		PlotWithBelts()
			:cCanvas(NULL),pt(NULL),gGreen(NULL),
			gYellow(NULL),lineOne(NULL),gr(NULL),
			legend(NULL),grobs(NULL),ptCMSPreli(NULL)
		{};
		~PlotWithBelts();
		PlotWithBelts(
				double *r1sigLt, double *r1sigHt, 
				double *r2sigLt, double *r2sigHt, 
				double *rmeant,  int npointst,  
				double *xpointst, string ssavetmp, TPaveText *pttmp,
				double xmintmp =100, double xmaxtmp =260, double ymintmp=0.5, double ymaxtmp=100, bool logYtmp=true,
				string stitlet ="; Higgs mass, m_{H} (GeV/c^{2}); Exclusion Limit, r = #sigma_{95%}/#sigma_{SM}");
		PlotWithBelts(
				double *r1sigLt, double *r1sigHt, 
				double *r2sigLt, double *r2sigHt, 
				double *rmeant,  double *robst, int npointst,  
				double *xpointst, string ssavetmp, TPaveText *pttmp,
				double xmintmp =100, double xmaxtmp =260, double ymintmp=0.5, double ymaxtmp=100, bool logYtmp=true,
				string stitlet ="; Higgs mass, m_{H} (GeV/c^{2}); Exclusion Limit, r = #sigma_{95%}/#sigma_{SM}");
		TGraph *getGreenGraph(){return gGreen;};
		TGraph *getYellowGraph(){return gYellow;};
		TGraph *getMeanGraph(){return gr;};
		TGraph *getObsGraph(){return grobs;};
		TLine *getLine(){return lineOne;};
		TPaveText *getTitlePT(){return pt;};
		TPaveText *getCMSPT(){return ptCMSPreli;};

		TCanvas *getCanvas(){return cCanvas;};
		TLegend *getLegend(){return legend;};
		void setSavePath(string tmp){ssave=tmp;};
		void drawLegend(string s_gr, string s_gGreen, string s_gYellow, string s_grobs="");
		void plot();
		void save();
		void clear();

	private:
		TCanvas *cCanvas;
		TLegend *legend;
		bool logY;
		string stitle;
		string ssave;

		double xmin, xmax, ymin, ymax;
		TPaveText *pt;
		TGraph *gGreen;
		TGraph *gYellow;
		TLine *lineOne;
		TGraph *gr;
		TGraph *grobs;
		TPaveText *ptCMSPreli;

		int bDrawGrobs;
		double *r1sigL, *r1sigH, *r2sigL, *r2sigH, *rmean, *robs, *xpoints;
		int npoints;

};
class DrawSigBkgPdfs{
	public:
		DrawSigBkgPdfs(double (*pdfs)(double *, double *), double *pars, int npars,
			double (*pdfb)(double *, double *), double *parb, int nparb, double Ns, double Nb,
				double xstart, double xstop, bool logY, string ssave, string stitle, bool debug=0);
		~DrawSigBkgPdfs();
		void draw();
		void save();
		TCanvas *getCanvas(){return cCanvas;};
		TLegend *getLegend(){return legend;};
		void setSavePath(string tmp){_ssave=tmp;};
		void drawLegend(string s_tot, string s_s, string s_b);
		void setLogY(bool b){_logY=b;};
		void setDebug(bool b){_debug=b;};

		TH1F *getHistS(){return hs;};
		TH1F *getHistB(){return hb;};
		TH1F *getHistTot(){return htot;};
		TF1 *getTf1S(){return fs;};
		TF1 *getTf1B(){return fb;};
	private:
		TCanvas *cCanvas;
		string _ssave;
		string _stitle;
		bool _logY;
		TLegend *legend;
		bool _debug;

		TF1 *fs;
		TF1 *fb;
		double (*_pdfs)(double *x, double *par); // normalized to unity
		double (*_pdfb)(double *x, double *par); // normalized to unity
		double *_pars;
		double *_parb;
		int _npars, _nparb;
		double _Ns, _Nb;		
		double _xstart, _xstop;
		TH1F *hs, *hb, *htot;
};
class DrawEvolution2D{
	public:
		DrawEvolution2D(vector<double> vx, vector<double> vy, string stitle, string ssave, TPaveText *pt, bool debug=0);
		~DrawEvolution2D();

		void setDrawPosteriorPdf(double limit){_drawPosteriorPdf = true; _limit = limit;};
		void draw();
		void save();
		TCanvas *getCanvas(){return cCanvas;};
		TLegend *getLegend(){return legend;};
		void setSavePath(string tmp){_ssave=tmp;};
		void drawLegend(string s_tot, string s_s, string s_b){};
		void setLogY(bool b){_logY=b;};
		void setDebug(bool b){_debug=b;};
		TPaveText *getTitlePT(){return _pt;};
		TGraph *getGraph(){return graph;};
		TLine *getLine(){return lineOne;};
	private:
		TCanvas *cCanvas;
		string _ssave;
		string _stitle;
		bool _logY;
		TLegend *legend;
		bool _debug;
		TPaveText *_pt;

		vector<double> _vx, _vy;
		
		bool _drawPosteriorPdf;
		double _limit;

		TLine *lineOne;
		TGraph *graph;
};
class DrawTGraph2D{
	public:
		DrawTGraph2D(vector< std::pair<double, double> > vxy, vector<double> vz, string stitle, string ssave, TPaveText *pt, bool debug=0);
		~DrawTGraph2D();

		void draw();
		void save();
		TCanvas *getCanvas(){return cCanvas;};
		TLegend *getLegend(){return legend;};
		void setSavePath(string tmp){_ssave=tmp;};
		void drawLegend(string s_tot, string s_s, string s_b){};
		void setLogY(bool b){_logY=b;};
		void setDebug(bool b){_debug=b;};
		TPaveText *getTitlePT(){return _pt;};
		TGraph2D *getGraph(){return graph;};
		TLine *getLine(){return lineOne;};
	private:
		TCanvas *cCanvas;
		string _ssave;
		string _stitle;
		bool _logY;
		TLegend *legend;
		bool _debug;
		TPaveText *_pt;

		vector<double> _vx, _vy, _vz;
		
		TLine *lineOne;
		TGraph2D *graph;
};

class PlotXvsCummulativeProb{
	public:
		PlotXvsCummulativeProb(
				double* vx_not_sorted, int nexps,
				string ssave, string stitle,TPaveText *pt);
		PlotXvsCummulativeProb(
				vector<double> vx_not_sorted,
				string ssave, string stitle,TPaveText *pt);
		PlotXvsCummulativeProb(
				vector<double> vx_sorted,
				vector<double> vcumulaP_sorted,
				double m1s, double p1s, double m2s, double p2s,
				string ssave, string stitle,TPaveText *pt);
		~PlotXvsCummulativeProb();
		void draw();
		void save();
		TCanvas *getCanvas(){return cCanvas;};
		TLegend *getLegend(){return legend;};
		void setSavePath(string tmp){_ssave=tmp;};
		void drawLegend(string s_tot, string s_s, string s_b){};
		void setLogY(bool b){_logY=b;};
		void setDebug(bool b){_debug=b;};
		TPaveText *getTitlePT(){return _pt;};
		TGraphErrors *getGraph(){return graph;};
		TLine *getLine(){return lineOne;};
		TH1F *getDifferencialHist(){return diffHist;};
		void setGraphDrawOption(string so){_graphOption=so;};

		// --- to show the error bar of bins
		void setTotalExps(int nexps){_nexps=nexps;};
		void setShowErrorBar(bool b){bShowErrorBar=b;}

		void drawDifferencial(string stitle, string ssave);

	private:
		TCanvas *cCanvas;
		string _ssave;
		string _stitle;
		bool _logY;
		TLegend *legend;
		bool _debug;
		TPaveText *_pt;

		vector<double> _vx, _vy;

		TLine *lineOne;
		TGraphErrors *graph;
		TH1F *diffHist;

		double _1SigmaLow, _1SigmaHigh, _2SigmaLow, _2SigmaHigh;

		bool bShowErrorBar;
		int _nexps;
	
		string _graphOption;	
};
class DrawPdfM2logQ{
	public:
		DrawPdfM2logQ(vector<double> vPdfM2logQ_sb, vector<double> vPdfM2logQ_b, 
				double m2lnQ_d, string sdata, string stitle, string ssave, TPaveText *pt); // show the famous -2lnQ with s, b, d
		DrawPdfM2logQ(vector< pair<double,double> > vPdfM2logQ_sb, vector< pair<double,double> > vPdfM2logQ_b, 
				double m2lnQ_d, string sdata, string stitle, string ssave, TPaveText *pt); // show the famous -2lnQ with s, b, d
		~DrawPdfM2logQ();
		void draw(); void draw2Hist(); void draw2Graph();
		TCanvas *getCanvas(){return cCanvas;};
		TLegend *getLegend(){return legend;};
		void setSavePath(string tmp){_ssave=tmp;};
		void drawLegend(string s_sb, string s_b, string s_d);
		void setLogY(bool b){_logY=b;};
		void setDebug(bool b){_debug=b;};
		TPaveText *getTitlePT(){return _pt;};
		TLine *getLine(){return lineOne;};
		TGraph *getSignalGraph(){return gSB;};
		TGraph *getBackgroundGraph(){return gB;};
		TH1F *getSignalHist(){return hSB;};
		TH1F *getBackgroundHist(){return hB;};
	private:
		TCanvas *cCanvas;
		string _ssave;
		string _stitle;
		bool _logY;
		TLegend *legend;
		bool _debug;
		TPaveText *_pt;
		TLine *lineOne;

		TGraph *gSB, *gB;
		TH1F *hSB, *hB;

		string _lineTitle; 
		double _m2lnQ_d;
		vector<double> _vsb, _vb; //inputs from pseudo experiments
		vector< pair<double,double> > _vpsb, _vpb; // inputs from analytic formula
		// tag the inputs are _vsb or _vpsb 
		bool bFromPseudoExps;
};

const int maxNumGraphs=10;
const int _glineStyle[maxNumGraphs]={2,2,8,8,1,1, 1,1,1,1};
const int _glineColor[maxNumGraphs]={4,2,6,1,5,3, 3, 3, 3, 3};
const int _gmarkerStyle[maxNumGraphs]={20,24,3,26,22,23, 23,23,23,23};
class DrawMultiGraph{
	public:
		DrawMultiGraph(double *x1, double *y1, int n1, string title1,
				string stitle, string ssave, TPaveText *pt,
				double xmin, double xmax, double ymin, double ymax, bool logY);
		DrawMultiGraph(double *x1, double *y1, int n1, string title1,
				double *x2, double *y2, int n2, string title2, string stitle, string ssave, TPaveText *pt,
				double xmin, double xmax, double ymin, double ymax, bool logY);
		~DrawMultiGraph();
		void add(double *x, double *y, int n, string title);
		void draw();
		TGraph* getGraph(int i){if(i<maxNumGraphs && i>=0) return _gr[i];else return 0;};

		TLegend *getLegend(){
			//cout<<"delete me2 legendY2="<<legend->GetY2()<<endl;
		return legend;};
		void setSavePath(string tmp){_ssave=tmp;};
		void drawLegend(string s_sb, string s_b, string s_d);
		void setLogY(bool b){_logY=b;};
		void setDebug(bool b){_debug=b;};
		TPaveText *getTitlePT(){return _pt;};
		TLine *getLine(){return lineOne;};
		TCanvas *getCanvas(){return cCanvas;};
		void drawOption(int i, string op){if(i<maxNumGraphs && i>=0)_drawOptions[i]=op;};
		void setLegendXY(double x, double y){legendX=x; legendY=y;};
	private:
		TCanvas *cCanvas;
		string _ssave;
		string _stitle;
		bool _logY;
		TLegend *legend;
		bool _debug;
		TPaveText *_pt;
		TLine *lineOne;

		TGraph *_gr[maxNumGraphs];//currenly limit to maxNumGraphs graphs..
		int _numGraphAdded;
		double *_x[maxNumGraphs], *_y[maxNumGraphs];
		int _n[maxNumGraphs];
		string _title[maxNumGraphs];
		double _xmin,_xmax,_ymin,_ymax;
		string _drawOptions[maxNumGraphs];
		double legendX, legendY;
};
class DrawPdfRLikelihood{
	public:
		DrawPdfRLikelihood(double r95, double *par, int npar, string stitle, string ssave, TPaveText *pt);
		~DrawPdfRLikelihood();
		void draw();
		TH1F* getHist(){return hist;};

		TLegend *getLegend(){return legend;};
		void setSavePath(string tmp){_ssave=tmp;};
		void drawLegend(string s_hist, string s_arrow){};
		void setLogY(bool b){_logY=b;};
		void setDebug(bool b){_debug=b;};
		TPaveText *getTitlePT(){return _pt;};
		TArrow *getArrow(){return arrow;};
		
	private:
		TCanvas *cCanvas;
		string _ssave;
		string _stitle;
		bool _logY;
		TLegend *legend;
		bool _debug;
		TPaveText *_pt;
		TArrow *arrow;
		TH1F *hist;
		double _r95;
		double *_par;
		int _npar;
};
};
#endif   /* ----- #ifndef PLOTUTILITIES_H  ----- */

