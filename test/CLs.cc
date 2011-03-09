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
using std::cout;
using std::endl;
using namespace lands;
void processParameters(int argc, const char* argv[]);
int debug=0; int nexps=100000; double s=0;  double b=0; double s_err = 0; double b_err = 0; double d=0;
int seed =1234; int pdftypeEs = 1; int pdftypeEb = 1; int EsEb_correlated = 0; int calcExpectedMeanLimit=0; int calcSignificance=0;	 
int main(int argc, const char* argv[]){
	processParameters(argc, argv);

	FloridaStyle();

	CRandom *rdm = new CRandom(seed);  //initilize a random generator

	// set up couting model
	CountingModel *cms=new CountingModel();
	cms->SetRdm(rdm);
	cms->AddChannel(s,b); // adding first channel
	cms->AddObservedData(0,d); //for the first channel
	// AddUncertainty(index_channel,  index_sample, err_in_relative_fraction,  pdf_type, index_correlation)
	//                 channel         signal        err            LogNormal               index_correlation
	cms->AddUncertainty(0,              0,            s_err,         pdftypeEs,                       1                ); 
	cms->AddUncertainty(0,              1,            b_err,         pdftypeEb,                       EsEb_correlated==0?2:1        ); 

	if(s_err !=0 || b_err !=0 ) { // decide to use or not systematic errors
		cms->SetUseSystematicErrors(true);
	}


	// initialize the calculator
	CLsBase frequentist;
	frequentist.SetDebug(debug);
	frequentist.SetRdm(rdm);

	frequentist.BuildM2lnQ(cms,nexps);
	double errs, errb, errsb;
	double cls = frequentist.CLs(errs);
	double clsb = frequentist.CLsb(errsb);
	double clb = frequentist.CLb(errb);
	cout<<"------------------------------------------------------------"<<endl;
	cout<<"Observed CLs = "<<cls<<" +/- "<<errs<<endl;
	cout<<"Observed CLsb = "<<clsb<<" +/- "<<errsb<<endl;
	cout<<"Observed CLb = "<<clb<<" +/- "<<errb<<endl;
	cout<<"------------------------------------------------------------"<<endl;

	CLsLimit clsr95;
	clsr95.SetDebug(debug);
	double rtmp;
	clsr95.SetAlpha(0.05);
	rtmp = clsr95.LimitOnSignalScaleFactor(cms, &frequentist,nexps);
	cout<<"------------------------------------------------------------"<<endl;
	cout<<"Observed Upper Limit on the ratio R at 95\% CL = "<<rtmp<<" +/- "<<clsr95.LimitErr()<<endl;
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

		TPaveText *pt = SetTPaveText(0.2, 0.7, 0.3, 0.9);
		pt = SetTPaveText(0.67, 0.4, 0.8, 0.7);
		pt->AddText("CLs statistical bands");
		string ssave="plot_cump_vs_r95";
		string stitle="; CLs Limit, r = #sigma_{95%}/#sigma_{SM}; cumulative probability;";
		PlotXvsCummulativeProb plotRvsP(all_calculated_R95s, cummulativeProbabilities,
				rm1s, rp1s, rm2s, rp2s,ssave, stitle, pt);
		plotRvsP.draw();
		

	}

	if(calcSignificance){
		int ntoysToDoSignificance = 10000000;
		cms->AddObservedData(0,d); //reset the observed data
		cout<<"\n Running "<<ntoysToDoSignificance<<" toys to evaluate significance for data "<<endl;
		double signi = frequentist.SignificanceForData(ntoysToDoSignificance);
		cout<<"------------------------------------------------------------"<<endl;
		cout<<" Observed Significance for the data = "<<signi<<endl;
		cout<<"------------------------------------------------------------"<<endl;

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
	int npar=5;
	if(argc>=6){
		s=atof( argv[1] );
		s_err=atof( argv[2] );
		b=atof( argv[3] );
		b_err=atof( argv[4] );
		d=atof( argv[5] );
		npar++;
		if(argc>=npar+1){
			calcExpectedMeanLimit=atoi(argv[npar]);
			npar++;
			if(argc>=npar+1){
				calcSignificance=atoi( argv[npar] );			
				npar++;
				if(argc>=npar+1){
					pdftypeEs=atoi( argv[npar] );			
					npar++;
					if(argc>=npar+1){
						pdftypeEb=atoi( argv[npar] );			
						npar++;
						if(argc>=npar+1){
							EsEb_correlated=atoi( argv[npar] );			
							npar++;
							if(argc>=npar+1){
								nexps=atoi( argv[npar] );			
								npar++;
								if(argc>=npar+1){
									seed=atoi( argv[npar] );			
									npar++;
									if(argc>=npar+1){
										debug=atoi( argv[npar] );			
									}
								}
							}
						}
					}
				}
			}
		}
	}
	if(argc<6) {
		cout<<"please use following format:"<<endl;
		cout<<"./CLs.exe signal err_sig_relative background err_bkg_relative data"<<endl;
		cout<<endl;
		cout<<" or "<<endl;
		cout<<"./CLs.exe signal err_sig_relative background err_bkg_relative data calcExpectedMeanLimit=0 calcSignificance=0 pdftypeEs=1 pdftypeEb=1 EsEb_correlated=0 ntoys=100000 seed=1234 debug=0"<<endl;
		cout<<endl<<" detail of the parameters: "<<endl;
		cout<<" signal: 		number of expected signal events " <<endl;
		cout<<" err_sig_relative: 	uncertainty on the signal, in relative fraction " <<endl;
		cout<<" background: 		number of expected background events " <<endl;
		cout<<" err_bkg_relative: 	uncertainty on the background, in relative fraction " <<endl;
		cout<<" data: 			number of observed events"<<endl;
		cout<<" calcExpectedMeanLimit:  calc the mean values of expected limit regardless of the data if 1 " <<endl;
		cout<<" calcSignificance:	calc the observed significance and expected mean value if 1 "<<endl;
		cout<<" pdftypeEs:		prior pdf on uncertainty of signal,  1 for LogNormal, 2 for TruncatedGaussian"<<endl;
		cout<<" pdftypeEb:		prior pdf on uncertainty of background,  1 for LogNormal, 2 for TruncatedGaussian"<<endl;
		cout<<" EsEb_correlated: 	100\% correlated between Es and Eb if 1,  uncorrelated if 0"<<endl;
		cout<<" ntoys: 			number of toy experiments to build -2lnQ "<<endl;
		cout<<" seed:			seed for random number generator "<<endl;
		cout<<" debug: 			debug level "<<endl;
		exit(0);
	}
	// print out  parameters  you just configured:
	cout<<"\n\t Your inputs are the following: "<<endl;	
	cout<<"**********************************************"<<endl;
	cout<<" expected signal yield                      = "<< s<<"+/-"<<s*s_err<<endl;
	cout<<" expected background yield                  = "<< b<<"+/-"<<b*b_err<<endl;
	cout<<" observed data events                       = "<< d<<endl;
	if(pdftypeEs!=1 && pdftypeEs!=2){ cout<<"pdftypeEs should be 1 or 2"<<endl; exit(0); }
	cout<<" prior pdf on uncertainty of signal is        "<<(pdftypeEs==1?"LogNormal":"TruncatedGaussian")<<endl; 
	if(pdftypeEb!=1 && pdftypeEb!=2){ cout<<"pdftypeEb should be 1 or 2"<<endl; exit(0); }
	cout<<" prior pdf on uncertainty of background is    "<<(pdftypeEb==1?"LogNormal":"TruncatedGaussian")<<endl; 
	cout<<" uncertainties of signal and background are "<< (EsEb_correlated==0?"uncorrelated":"100\% correlated")<<endl;
	cout<<(calcExpectedMeanLimit==0?" NOT":"")<<" do calculation of mean limit"<<endl;
	cout<<(calcSignificance==0?" Not":"")<<" do significance calculation"<<endl;
	cout<<" seed = "<< seed<<endl;
	cout<<" debuglevel = "<<debug<<endl;
	cout<<"**********************************************"<<endl;
}
