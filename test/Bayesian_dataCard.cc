#include <iostream>
//#include "CLsLimit.h"
#include "CRandom.h"
#include "PlotUtilities.h"
#include "FloridaStyle.C"
#include <TError.h>
#include "TSystem.h"
#include "CountingModel.h"
#include "BayesianBase.h"
#include "LimitBands.h"
#include "UtilsROOT.h"
using std::cout;
using std::endl;
using namespace lands;
void processParameters(int argc, const char* argv[]);
int debug=0; int nexps=20000; double s=0;  double b=0; double s_err = 0; double b_err = 0; int d=0;
int seed =1234; int pdftypeEs = 1; int pdftypeEb = 1; int EsEb_correlated = 0; int calcExpectedMeanLimit=0; //int calcSignificance=0;	 
int XSprior = 10;
int asimov = -1; // -1 for using the input dataset whatever you provide,  0 for using asimov bkg only dataset, 1 for asimov sig+bkg
const char * fileName;
int main(int argc, const char* argv[]){
	processParameters(argc, argv);

	FloridaStyle();

	CRandom *rdm = new CRandom(seed);  //initilize a random generator

	// set up couting model
	CountingModel *cms=new CountingModel();
	cms->SetRdm(rdm);

	ConfigureModel(cms, fileName); 
	cms->SetUseSystematicErrors(true);

	cms->UseAsimovData(asimov);

	cms->Print(debug);
	cms->RemoveChannelsWithExpectedSignal0orBkg0();

	cms->Print(debug);
	
	BayesianBase bys(cms, 0.05, 1.e-3);
	bys.SetNumToys(nexps);
	bys.SetDebug(debug);
	if(XSprior==20)bys.SetCrossSectionPrior(corr);
	else bys.SetCrossSectionPrior(flat);
	double rtmp;
	rtmp = bys.Limit();
	cout<<"------------------------------------------------------------"<<endl;
	cout<<"Observed Upper Limit on the ratio R at 95\% CL = "<<rtmp<<endl;
	cout<<"------------------------------------------------------------"<<endl;

	// draw 

	/*
	bys.PosteriorPdf();
	TPaveText *pt = SetTPaveText(0.2, 0.7, 0.3, 0.9);
	//start_time=cur_time; cur_time=clock(); cout << "\t TIME_in_BayesianLimit ppdf: "<< (cur_time - start_time)/1000. << " millisec\n";
	DrawEvolution2D pdfr(bys.GetVR(), bys.GetVP(), ";r;Likelihood", "likelihood_pdf", pt);	
	pdfr.draw();
	pdfr.getLine()->Delete();
	pdfr.getGraph()->SetMarkerSize(0.01);
	pdfr.save();
	*/
	

	cms->SetSignalScaleFactor(1.);
	if(calcExpectedMeanLimit){
		LimitBands lb(&bys, cms);	
		lb.SetDebug(debug);
		int noutcomes = 1000;
		lb.BysLimitBands(0.05, noutcomes, nexps);

		double rmean, rm1s, rm2s, rp1s, rp2s;	
		rmean=lb.GetBysLimitMean();
		rm1s=lb.GetBysLimit(-1);rm2s=lb.GetBysLimit(-2);rp1s=lb.GetBysLimit(1);rp2s=lb.GetBysLimit(2);

		// plot the distribution of limits of noutcomes
		// and the cummulative pdf
		TString ts(fileName); ts+="_byslimits";
		FillTree(ts, lb.GetDifferentialLimitsBys());

		cout<<"------------------------------------------------------------"<<endl;
		cout<<"Expected R@95%CL (from -2sigma -1sigma  mean  +1sigma  +2sigma) : "<<endl;
		printf("--------- %10.5f %10.5f %10.5f %10.5f %10.5f\n", rm2s, rm1s, rmean, rp1s, rp2s);
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
							XSprior=atoi( argv[npar] );			
							if(XSprior !=10 and XSprior!=20) {
								cout<<"prior on the cross section should be \"10\" or \"20\""<<endl;
								exit(0);
							}
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
	if(argc<2) {
		cout<<"please use following format:"<<endl;
		cout<<"./Bayesian.exe inputFileName"<<endl;
		cout<<endl;
		cout<<" or "<<endl;
		cout<<"./Bayesian.exe inputFileName calcExpectedMeanLimit=0 ntoys=100000 seed=1234 debug=0 XSprior=10"<<endl;
		cout<<endl<<" detail of the parameters: "<<endl;
		cout<<" inputFileName:          The input data card designed by Andrey Korytov "<<endl;
		cout<<" calcExpectedMeanLimit:  calc the mean values of expected limit regardless of the data if 1 " <<endl;
		cout<<" ntoys: 			number of toy experiments to avarage out uncertainties "<<endl;
		cout<<" seed:			seed for random number generator "<<endl;
		cout<<" debug: 			debug level "<<endl;
		cout<<" XSprior:                prior on the cross section ratio,  10 for flat prior, 20 for reference prior " <<endl;
		cout<<" asimov:                 0 for bkg only, 1 for sig+bkg, others for the input whatever you provide "<<endl;
		exit(0);
	}
	// print out  parameters  you just configured:
	cout<<"\n\t Your inputs are the following: "<<endl;	
	cout<<"**********************************************"<<endl;
	cout<<" Input data card is:   "<<fileName<<endl; 
	cout<<" prior on the cross section ratio = "<<(XSprior==10?"flat":"reference")<<endl;
	cout<<(calcExpectedMeanLimit==0?" NOT":"")<<" do calculation of mean limit"<<endl;
	cout<<" seed = "<< seed<<endl;
	cout<<" debuglevel = "<<debug<<endl;

	if(asimov==0) cout<<" Using Asimov bkg only dataset"<<endl;
	if(asimov==1) cout<<" Using Asimov signal + bkg  dataset"<<endl;

	cout<<"**********************************************"<<endl;
}
