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
int debug=0; int nexps=1000; 
int seed =1234; int calcExpectedMeanLimit=0; 
int asimov = -1; // -1 for using the input dataset whatever you provide,  0 for using asimov bkg only dataset, 1 for asimov sig+bkg
int oneside = 1; 
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
	cms_global= cms;
	vdata_global=cms->Get_v_data();

	double r95;
	double tmp;
/*
	//double x0 =  MinuitFit(3, tmp, tmp);
	double x0 =  MinuitFit(2, tmp, tmp);
	cout<<"x0 = "<<x0<<endl;
	if(myMinuit) cout<<"myMinuit != 0 "<<endl;
	x0 =  MinuitFit(3, tmp, tmp);
	myMinuit->Release(0);
	myMinuit->Command("SCAn 1 10 0 0.3");
	TCanvas c("ct","ct");
	TGraph *gr = (TGraph*)myMinuit->GetPlot();
	//TH1F *h = (TH1F*)gr->GetHistogram();
	//h->Scale(0.5);
	//h->SetMarkerStyle(21);
	//h->Draw("lp");
	gr->SetMarkerStyle(21);
	gr->Draw("alp");
	c.SaveAs("del.root");
	c.SaveAs("del.gif");
	gr->Print("");

	double x1 =  MinuitFit(3, tmp, tmp, 1);
	cout<<"x1 = "<<x1<<endl;
	if(myMinuit) cout<<"myMinuit != 0 "<<endl;
	myMinuit->Release(0);
	myMinuit->Command("SCAn 1 10 0 0.3");
	TCanvas c1("ct1","ct1");
	TGraph *gr1 = (TGraph*)myMinuit->GetPlot();
	gr1->SetMarkerStyle(21);
	gr1->Draw("alp");
	c1.SaveAs("del1.root");
	c1.SaveAs("del1.gif");

	gr1->Print("");
  */  
	double tmperr;
	double y0_1 =  MinuitFit(3, tmp, tmperr, 0);
	cout<<y0_1<<" fitter u="<<tmp<<" +/- "<<tmperr<<endl;
	double tmpr = 0;
	double y0_2 =  MinuitFit(2, tmpr, tmperr) ;
	cout<<y0_2<<" fitter u="<<tmp<<" +/- "<<tmperr<<endl;
	
	double x1 =0, x2 =1;
	double y0 = y0_1;  
	if ( (y0_1 > y0_2) && tmpr>0) {
		y0 = y0_2;
		x1 = tmpr; x2=2*tmpr;
	}


	double y =  MinuitFit(3, tmp, tmp, x2 );

	//------------
	//If we have a background as a Gaussian distribution G(x|b) with mean b,
	//sigma sqrt(b), and observe x = b, then the upper limit on signal is
	//mu = 1.64*sqrt(b), which gives a 5% chance for P(obs < b).
	//
	//Therefore:
	//
	//1) muhat = 0, given the observation x=b
	//
	//2) -2 * ln (lambda(mu)) = -2 * ln( G(x|b+mu) / G(x|b+muhat) ) = 1.64^2
	//
	//If we use the 1.96-rule, it may *artificially* improve, but not cure,
	//the coverage for low statistics case, but would now give wrong coverage
	//in asymptotic.

	//http://en.wikipedia.org/wiki/Chi-square_distribution
	//double CI = 1.921;  // = 1.96**2/2 ,  probably for two sided 
	//double CI = 1.64*1.64/2.;
	double CI ;
	if (oneside==0) CI= 1.921;  // = 1.96**2/2 ,  probably for two sided 
	else CI = 1.64*1.64/2.;
	double precision = 0.001;
	int nsearched = 2;
	//1.925 for 95% CL,   ....   
	//              assume the profile likelihood function to be an increasing function
	//              bisection search ...
	//              y1 always > y0

	if(fabs((y-y0)/2. - CI) > precision ){
		//first, got a number > 1.921,  otherwise increase it by a factor of 10 ...
		while( (y-y0)/2. < CI ){
			x1 =  x2;
			x2 *=10.;
			y =  MinuitFit(3, tmp, tmp, x2 );
			nsearched++;
		}
		y = MinuitFit(3, tmp, tmp, (x1+x2)/2. );
		while( fabs((y-y0)/2. - CI)/CI > precision ){
			double  tmpx = (x1+x2)/2.;
			if( (y-y0)/2. < CI ){
				x1 = tmpx; 
			}else{
				x2 = tmpx;
			}	
			y = MinuitFit(3, tmp, tmp, (x1+x2)/2. );
			nsearched++;
		}
		r95 = (x2+x1)/2.;
	}else{
		r95 = x2;
	}

	cout<<"r95 = "<<r95<<",  "<<nsearched<<" steps"<<endl;



	if(debug>=10){ // show a plot for  -log(Lambda(mu)) vs. mu 
		//double x0 =  MinuitFit(3, tmp, tmp, 0 );
		double x0 = y0;
		cout<<"x0 = "<<x0<<endl;
		vector<double> vrxsec, vmlnq; 
		vrxsec.clear(); vmlnq.clear();
		for(double r=0; r<=2*r95; r+=r95/10.){
			vrxsec.push_back(r);
			double x3 = MinuitFit(3,tmp, tmp, r) ;
			cout<<"r= "<<r<<", x3 = "<<x3<<endl;
			double m2lnQ = x3 - x0;
			cout<<"Profiled:  r ="<<r<<"	m2lnQ="<<m2lnQ<<endl;
			vmlnq.push_back(m2lnQ/2.0);
			if(m2lnQ/2.>3) break;
		}
		DrawEvolution2D profiledL(vrxsec, vmlnq, ";r = #sigma/#sigma_{SM}; -log#lambda(r)", "profiledL", pt);
		profiledL.setLogY(0);
		profiledL.draw();
	}

}
void processParameters(int argc, const char* argv[]){
	int npar=1;
	if(argc>=2){
		fileName = argv[1];
		npar++;
		if(argc>=npar+1){
			oneside=atoi( argv[npar] );			
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
		cout<<"./ProfileLikelihoodApproxLimit.exe inputFileName"<<endl;
		cout<<endl;
		cout<<" or "<<endl;
		cout<<"./ProfileLikelihoodApproxLimit.exe inputFileName oneside=1 calcExpectedMeanLimit=0 ntoys=1000 seed=1234 debug=0 asimov=-1"<<endl;
		cout<<endl<<" detail of the parameters: "<<endl;
		cout<<" inputFileName:          The input data card designed by Andrey Korytov "<<endl;
		cout<<" oneside:                0 for two-sided upper limit,  1 (default) for one-sided "<<endl;
		cout<<" calcExpectedMeanLimit:  calc the mean values of expected limit regardless of the data if 1 " <<endl;
		cout<<" ntoys: 			number of toy experiments to build bands"<<endl;
		cout<<" seed:			seed for random number generator "<<endl;
		cout<<" debug: 			debug level "<<endl;
		cout<<" asimov:                 0 for bkg only, 1 for sig+bkg, others for the input whatever you provide "<<endl;
		exit(0);
	}
	// print out  parameters  you just configured:
	cout<<"\n\t Your inputs are the following: "<<endl;	
	cout<<"**********************************************"<<endl;
	cout<<" Input data card is:   "<<fileName<<endl; 
	cout<<(oneside==0?" two":"one")<<"-sided upper limit"<<endl;
	cout<<(calcExpectedMeanLimit==0?" NOT":"")<<" do calculation of mean limit"<<endl;
	cout<<" seed = "<< seed<<endl;
	cout<<" debuglevel = "<<debug<<endl;

	if(asimov==0) cout<<" Using Asimov bkg only dataset"<<endl;
	if(asimov==1) cout<<" Using Asimov signal + bkg  dataset"<<endl;

	cout<<"**********************************************"<<endl;
}
