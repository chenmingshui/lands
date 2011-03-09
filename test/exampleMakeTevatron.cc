#include <vector>
#include <iostream>
#include <iomanip>
#include <ctime>
#include <math.h>

#include <TH1F.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TError.h>
#include <TFile.h>
#include <TTree.h>
#include <TPaveText.h>

#include "FloridaStyle.C"
#include "PlotUtilities.h"


#include "TTree.h"
#include "TSystem.h"

using std::cout;
using std::endl;
using namespace stats;
char ctmp[255];
string ssave;



int debug = 0;
int main(int argc, const char* argv[]){
	const int nMassPoints = 2;
	double mass_xaxis[nMassPoints] = {130, 140};
	double meanValues_yaxis[nMassPoints] = {1.2, 2.3};
	double positiveOneSigma[nMassPoints] = {1.28, 2.5};
	double negativeOneSigma[nMassPoints] = {1.10, 2.0};
	double positiveTwoSigma[nMassPoints] = {1.4, 3.0};
	double negativeTwoSigma[nMassPoints] = {0.95, 1.8};

	TString title = "";
	TString title_y = "sigma_{95%}/sigma_{SM}";
	TString title_x = "Higgs mass (GeV/c^{2})";
	double xaxis_min = 100;
	double xaxis_max = 200;
	double yaxis_min = 0.1;
	double yaxis_max = 5;
	bool yaxis_log = false;

	TString filename_to_save = "del.gif"

	if(!debug)gErrorIgnoreLevel=5000; // do not show the message when saving plots
	FloridaStyle(); // plots style


	TPaveText *pt;
	PlotWithBelts tev(
			negativeOneSigma, positiveOneSigma, negativeTwoSigma, positiveTwoSigma,
			meanValues_yaxis, meanValues_yaxis, nMassPoints,  
			mass_xaxis, filename_to_save.Data(), pt,
			xaxis_min, xaxis_max, yaxis_min, yaxis_max, yaxis_log, (title+title_x+title_y).Data());
	tev.plot();
	tev.getObsGraph()->SetMarkerStyle(1);//remove marker
	tev.drawLegend("CLs 95% CL exclusion: mean","95% CL exclusion: 68% band", "95% CL exclusion: 95% band", "95% CL exclusion: mean (no sys)");
	MoveLegend(tev.getLegend(),0.5,0.6);
	tev.getMeanGraph()->SetLineStyle(1);
	tev.getMeanGraph()->SetLineColor(kBlue);
	tev.getMeanGraph()->SetLineWidth(2);
	tev.getObsGraph()->SetLineStyle(2);
	tev.getObsGraph()->SetLineWidth(1);
	tev.getLine()->SetLineColor(kRed);
	tev.save();	

	return 1;
}
