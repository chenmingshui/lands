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
int debug=0; int nexps=1000000; 
int seed =1234; int calcSignificance=0;	 
int calcProfiledLikelihoodSigificance = 1;
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

	double tmp, tmperr;
	if(calcProfiledLikelihoodSigificance >=1){
		vdata_global=cms->Get_v_data();

		double x2 =  MinuitFit(2, tmp, tmperr);
		cout<<"fitted r = "<<tmp<<endl;
		double m2lnQ = MinuitFit(3,tmp, tmp) - x2;
		double sig_data = sqrt(fabs(m2lnQ));
		cout<<"Observed significance using PLR method = "<<sig_data<<endl;

		if(debug>=10) { // show a plot for   -log(Lambda(mu)) vs. mu ...
			for(double r=0; r<2; r+=0.1){
				m2lnQ = MinuitFit(3,tmp, tmp, r) -x2; 
			}
		}
	}


	if(calcProfiledLikelihoodSigificance >= 2){
		SignificanceBands lb(&frequentist, cms);	
		lb.SetDebug(debug);
		int noutcomes = 1000;
		lb.Bands(noutcomes);

		double rmean, rm1s, rm2s, rp1s, rp2s;	
		rmean=lb.GetSignificanceMean();
		rm1s=lb.GetSignificance(-1);rm2s=lb.GetSignificance(-2);rp1s=lb.GetSignificance(1);rp2s=lb.GetSignificance(2);
		cout<<"------------------------------------------------------------"<<endl;
		cout<<"Expected significance (from -2sigma -1sigma  mean  +1sigma  +2sigma) : "<<endl;
		printf("--------- %10.5f %10.5f %10.5f %10.5f %10.5f\n", rm2s, rm1s, rmean, rp1s, rp2s);
		cout<<"------------------------------------------------------------"<<endl;

		TString ts(fileName); ts+="_PLRsignificances";
		FillTree(ts, lb.GetDifferentialSignificances());

		vector<double> difflimits=lb.GetDifferentialSignificances();
		TCanvas *c=new TCanvas("csig","cSig");
		c->SetLogy(1);
		TH1F *h=new TH1F("h",";ProfiledLikelihood significance; entries", 200, 0, 8);	
		for(int i=0; i<difflimits.size(); i++) h->Fill(difflimits[i]);
		h->Draw();
		Save(c, "differential_PL_sig");

		vector<double> all_calculated_R95s;
		vector<double> cummulativeProbabilities; // all_calculated_R95s has been sorted, each of them has a cummulative probability
		SortAndCumulative(difflimits, all_calculated_R95s, cummulativeProbabilities);

		pt = SetTPaveText(0.67, 0.4, 0.8, 0.7);
		pt->AddText("PL significance statistical bands");
		string ssave="plot_cump_vs_PLsig";
		string stitle="; ProfiledLikelihood significance; cumulative probability;";
		PlotXvsCummulativeProb plotRvsP(all_calculated_R95s, cummulativeProbabilities,
				rm1s, rp1s, rm2s, rp2s,ssave, stitle, pt);
		plotRvsP.draw();


	}


	cms->SetSignalScaleFactor(1.);
	int ntoysToDoSignificance = nexps; //10000000;
	frequentist.SetModel(cms);
	if(calcSignificance==1){
		//cms->AddObservedData(0,d); //reset the observed data
		cout<<"\n Running "<<ntoysToDoSignificance<<" toys to evaluate significance for data "<<endl;
		double signi = frequentist.SignificanceForData(ntoysToDoSignificance);
		if(debug){
			cout<<"Q_b_data = "<<frequentist.Get_m2lnQ_data()<<endl;	
			vector<double> vclb = frequentist.Get_m2logQ_b();
			TString  s = fileName; 
			s+="_hybridSig_ts"; s+=testStatistics;
			s+="_seed"; s+=seed;
			FillTree(s, vclb);
			
		}
		cout<<"------------------------------------------------------------"<<endl;
		cout<<" Observed Significance for the data = "<<signi<<endl;
		cout<<"------------------------------------------------------------"<<endl;
	}
	if(calcSignificance>=2){
		vector<double> vsignificance;
		vector<double> vsignificance_cp;
		double significance[5];
		cout<<"\n Running "<<ntoysToDoSignificance<<" toys to evaluate significance for expected values "<<endl;
		double significance_mean = frequentist.SignificanceComputation(10000,  ntoysToDoSignificance, vsignificance, vsignificance_cp );
		GetBandsByLinearInterpolation(vsignificance,vsignificance_cp, significance[1], significance[3], significance[0], significance[4] );
		significance[2]=GetBandByLinearInterpolation(vsignificance, vsignificance_cp, 0.5);
		cout<<"------------------------------------------------------------"<<endl;
		cout<<"Expected Significance (from -2sigma -1sigma  median  +1sigma  +2sigma  mean) : "<<endl;
		printf("---------- %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n",
				significance[0], significance[1], significance[2], significance[3], significance[4], significance_mean); 
		cout<<"------------------------------------------------------------"<<endl;
	}

	return 1;
}
void processParameters(int argc, const char* argv[]){
	int npar=1;
	if(argc>=2){
		fileName = argv[1];
		npar++;
		if(argc>=npar+1){
			calcProfiledLikelihoodSigificance=atoi(argv[npar]);
			npar++;
			if(argc>=npar+1){
				calcSignificance=atoi( argv[npar] );			
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
	}
	if(argc<2) {
		cout<<"please use following format:"<<endl;
		cout<<"./Significance_dataCard.exe inputFileName"<<endl;
		cout<<endl;
		cout<<" or "<<endl;
		cout<<"./Significance_dataCard.exe inputFileName calcPLRsignificance=1 calcFRQsignificance=0 ntoys=1000000 seed=1234 debug=0 testStatistics=1 rule=1"<<endl;
		cout<<endl<<" detail of the parameters: "<<endl;
		cout<<" inputFileName:          The input data card designed by Andrey Korytov "<<endl;
		cout<<" calcPLRsignificance:    1 calc the significance with PLR method for the observation, 2 calc expected mean value regardless of the data  " <<endl;
		cout<<" calcFRQsignificance:    1 calc the significance with frequentist method for the observation, 2 calc expected mean value regardless of the data  " <<endl;
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
	if(calcProfiledLikelihoodSigificance==1) cout<<"calc PLR significance for data"<<endl;
	if(calcProfiledLikelihoodSigificance==2) cout<<"calc PLR expected mean significance regardless of data"<<endl;
	if(calcSignificance==1) cout<<"calc Frequentist significance for data"<<endl;
	if(calcSignificance==2) cout<<"calc Frequentist expected mean significance regardless of data"<<endl;

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
