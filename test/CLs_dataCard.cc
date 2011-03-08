#include <iostream>
//#include <python2.4/Python.h>
#include "CLsLimit.h"
#include "CRandom.h"
#include "PlotUtilities.h"
#include "FloridaStyle.C"
#include <TError.h>
#include "TSystem.h"
#include "CountingModel.h"
#include "BayesianBase.h"
#include "LimitBands.h"
#include "SignificanceBands.h"
#include "UtilsROOT.h"
using std::cout;
using std::endl;
using namespace lands;
void processParameters(int argc, const char* argv[]);
int debug=0; int nexps=100000; double s=0;  double b=0; double s_err = 0; double b_err = 0; int d=0;
int seed =1234; int pdftypeEs = 1; int pdftypeEb = 1; int EsEb_correlated = 0; int calcExpectedMeanLimit=0; 
int testStatistics = 1, rule = 1; // default is CLs
int asimov = -1; // -1 for using the input dataset whatever you provide,  0 for using asimov bkg only dataset, 1 for asimov sig+bkg
const char* fileName;
int main(int argc, const char* argv[]){
	processParameters(argc, argv);

	FloridaStyle();
	TPaveText *pt = SetTPaveText(0.2, 0.7, 0.3, 0.9);

	CRandom *rdm = new CRandom(seed);  //initilize a random generator

	// set up couting model
	CountingModel *cms=new CountingModel();
	cms->SetRdm(rdm);
	ConfigureModel(cms, fileName); 
	cms->SetUseSystematicErrors(true);

	cms->UseAsimovData(asimov);
	cms->RemoveChannelsWithExpectedSignal0orBkg0();
	cms->Print();

	// initialize the calculator
	CLsBase frequentist;
	frequentist.SetDebug(debug);
	frequentist.SetRdm(rdm);
	frequentist.SetTestStatistics(testStatistics);
	cms_global= cms;
	//vdata_global=cms->Get_v_data();

	double tmp;
	vdata_global=cms->Get_v_data();

/*	double m2lnQ = MinuitFit(3,tmp, tmp) - MinuitFit(2, tmp, tmp);
	double sig_data = sqrt(fabs(m2lnQ));
	cout<<"Observed significance = "<<sig_data<<endl;
*/

	if(0){
		vector<double> vrxsec, vmlnq; 
		vrxsec.clear(); vmlnq.clear();
		for(double r=0; r<2; r+=0.01){
			vrxsec.push_back(r);
			double m2lnQ = MinuitFit(3,tmp, tmp, r) - MinuitFit(2, tmp, tmp);
			cout<<"Profiled:  r ="<<r<<"	m2lnQ="<<m2lnQ<<endl;
			vmlnq.push_back(m2lnQ/2.0);
		}
		DrawEvolution2D profiledL(vrxsec, vmlnq, ";rXsec; -log#lambda(rXsec)", "profiledL", pt);
		profiledL.draw();
	}

	frequentist.BuildM2lnQ(cms,nexps);
	double cls = frequentist.CLs();
	cout<<"------------------------------------------------------------"<<endl;
	cout<<"Observed CLs = "<<cls<<endl;
	cout<<"------------------------------------------------------------"<<endl;

	DrawPdfM2logQ pdfM2logQ(frequentist.Get_m2logQ_sb(),frequentist.Get_m2logQ_b(), frequentist.Get_m2lnQ_data(), 
			"-2lnQ on data", "; -2lnQ; entries", "lnq", pt);
	pdfM2logQ.draw();

	if(debug>=10) {
		cout<<"-2lnQ on data = "<<frequentist.Get_m2lnQ_data()<<endl;
		FillTree("m2lnQ_b.root", frequentist.Get_m2logQ_b());
		FillTree("m2lnQ_sb.root", frequentist.Get_m2logQ_sb());
		vector<double> qsb = frequentist.Get_m2logQ_sb();
		vector<double> qb = frequentist.Get_m2logQ_b();
		cout<<"-2lnQ for SB"<<endl;
		for(int i=0; i<qsb.size(); i++) cout<<qsb[i]<<endl;
		cout<<"-2lnQ for B"<<endl;
		for(int i=0; i<qb.size(); i++) cout<<qb[i]<<endl;
	}



	CLsLimit clsr95;
	clsr95.SetDebug(debug);
	clsr95.SetRule(rule);
	double rtmp;
	clsr95.SetAlpha(0.05);
	rtmp = clsr95.LimitOnSignalScaleFactor(cms, &frequentist,nexps);
	cout<<"------------------------------------------------------------"<<endl;
	cout<<"Observed Upper Limit on the ratio R at 95\% CL = "<<rtmp<<endl;
	cout<<"------------------------------------------------------------"<<endl;

	cms->SetSignalScaleFactor(1.);
	if(calcExpectedMeanLimit){
		LimitBands lb(&clsr95, &frequentist, cms);	
		lb.SetDebug(debug);
		int noutcomes = 1000;
		lb.CLsLimitBands(0.05, noutcomes, nexps);

		double rmean, rm1s, rm2s, rp1s, rp2s;	
		rmean=lb.GetCLsLimitMean();
		rm1s=lb.GetCLsLimit(-1);rm2s=lb.GetCLsLimit(-2);rp1s=lb.GetCLsLimit(1);rp2s=lb.GetCLsLimit(2);
		cout<<"------------------------------------------------------------"<<endl;
		cout<<"Expected R@95%CL (from -2sigma -1sigma  mean  +1sigma  +2sigma) : "<<endl;
		printf("--------- %10.5f %10.5f %10.5f %10.5f %10.5f\n", rm2s, rm1s, rmean, rp1s, rp2s);
		cout<<"------------------------------------------------------------"<<endl;

		TString ts(fileName); ts+="_clslimits";
		FillTree(ts, lb.GetDifferentialLimitsBys());

		vector<double> difflimits=lb.GetDifferentialLimitsCLs();
		TCanvas *c=new TCanvas("cme","cme");
		c->SetLogy(1);
		TH1F *h=new TH1F("h",";r=#frac{#sigma_{95%CL}}{#sigma_{SM}}; entries", 200, 0, 15);	
		for(int i=0; i<difflimits.size(); i++) h->Fill(difflimits[i]);
		h->Draw();
		Save(c, "differential_limits");

		vector<double> all_calculated_R95s;
		vector<double> cummulativeProbabilities; // all_calculated_R95s has been sorted, each of them has a cummulative probability
		SortAndCumulative(difflimits, all_calculated_R95s, cummulativeProbabilities);

		pt = SetTPaveText(0.67, 0.4, 0.8, 0.7);
		pt->AddText("CLs statistical bands");
		string ssave="plot_cump_vs_r95";
		string stitle="; CLs Limit, r = #sigma_{95%}/#sigma_{SM}; cumulative probability;";
		PlotXvsCummulativeProb plotRvsP(all_calculated_R95s, cummulativeProbabilities,
				rm1s, rp1s, rm2s, rp2s,ssave, stitle, pt);
		plotRvsP.draw();


	}

	return 1;
}
void processParameters(int argc, const char* argv[]){
	int npar=1;
	if(argc>=2){
		fileName = argv[1];
		npar++;
		if(argc>=npar+1){
			calcExpectedMeanLimit=atoi(argv[npar]);
			npar++;
			if(argc>=npar+1){
				nexps=atoi( argv[npar] );			
				npar++;
				if(argc>=npar+1){
					seed=atoi( argv[npar] );			
					npar++;
					if(argc>=npar+1){
						debug=atoi( argv[npar] );			
						npar++;
						if(argc>=npar+1){
							testStatistics=atoi( argv[npar] );			
							npar++;
							if(argc>=npar+1){
								rule=atoi( argv[npar] );			
								npar++;
								if(argc>=npar+1){
									asimov=atoi( argv[npar] );			
								}
							}
						}
					}
				}
			}
		}
	}
	if(argc<2) {
		cout<<"please use following format:"<<endl;
		cout<<"./CLs.exe inputFileName"<<endl;
		cout<<endl;
		cout<<" or "<<endl;
		cout<<"./CLs.exe inputFileName calcExpectedMeanLimit=0 ntoys=100000 seed=1234 debug=0 testStatistics=1 rule=1"<<endl;
		cout<<endl<<" detail of the parameters: "<<endl;
		cout<<" inputFileName:          The input data card designed by Andrey Korytov "<<endl;
		cout<<" calcExpectedMeanLimit:  calc the mean values of expected limit regardless of the data if 1 " <<endl;
		cout<<" ntoys: 			number of toy experiments to build -2lnQ "<<endl;
		cout<<" seed:			seed for random number generator "<<endl;
		cout<<" debug: 			debug level "<<endl;
		cout<<" testStatistics:		1 for Q_LEP, 2 for Q_TEV, 3 for Q_ATLAS "<<endl;
		cout<<" rule:                   1 for CLs,  2 for CLsb "<<endl;
		cout<<" asimov:                 0 for bkg only, 1 for sig+bkg, others for the input whatever you provide "<<endl;
		exit(0);
	}
	// print out  parameters  you just configured:
	cout<<"\n\t Your inputs are the following: "<<endl;	
	cout<<"**********************************************"<<endl;
	cout<<" Input data card is:   "<<fileName<<endl; 
	cout<<(calcExpectedMeanLimit==0?" NOT":"")<<" do calculation of mean limit"<<endl;
	cout<<" seed = "<< seed<<endl;
	cout<<" debuglevel = "<<debug<<endl;

	cout<<" testStatistics = "<<testStatistics<<", is ";
	if(testStatistics==2) cout<<" Tevatron type";
	else if(testStatistics==3) cout<<" ATLAS type,  Lamda(mu)";
	else cout<<" LEP type";
	cout<<endl;

	cout<<" rule = "<< rule<<", is ";
	if(rule==2) cout<<" CLsb";
	else cout<<" CLs";
	cout<<endl;

	if(asimov==0) cout<<" Using Asimov bkg only dataset"<<endl;
	if(asimov==1) cout<<" Using Asimov signal + bkg  dataset"<<endl;

	cout<<"**********************************************"<<endl;
}
