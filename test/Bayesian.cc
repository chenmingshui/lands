#include <iostream>
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
void processParameters(int argc, const char* argv[]);
int debug=0; int nexps=20000; double s=0;  double b=0; double s_err = 0; double b_err = 0; double d=0;
int seed =1234; int pdftypeEs = 1; int pdftypeEb = 1; int EsEb_correlated = 0; int calcExpectedMeanLimit=0; //int calcSignificance=0;	 
int prior = 10;
int main(int argc, const char* argv[]){
	processParameters(argc, argv);

	CRandom *rdm = new CRandom(seed);  //initilize a random generator

	// set up couting model
	CountingModel *cms=new CountingModel();
	cms->SetRdm(rdm);
	cms->AddChannel(s,b); // adding first channel
	cms->AddObservedData(0,d); //for the first channel
	// AddUncertainty(index_channel,  index_sample, err_in_relative_fraction,  pdf_type, index_correlation)
	//                 channel         signal        err            LogNormal               index_correlation
	cms->AddUncertainty(0,              0,            s_err,         pdftypeEs,                       "1"                ); 
	cms->AddUncertainty(0,              1,            b_err,         pdftypeEb,                       EsEb_correlated==0?"2":"1"        ); 

	if(s_err !=0 || b_err !=0 ) { // decide to use or not systematic errors
		cms->SetUseSystematicErrors(true);
	}



	BayesianBase bys(cms, 0.05, 1.e-3);
	
	bys.SetNumToys(nexps);
	bys.SetDebug(debug);
	cout<<"prior = "<<prior<<endl;
	if(prior==10)bys.SetCrossSectionPrior(flat);
	if(prior==20)bys.SetCrossSectionPrior(corr);
	if(prior==30)bys.SetCrossSectionPrior(prior_1overSqrtS);
	double rtmp;
	rtmp = bys.Limit();
	cout<<"------------------------------------------------------------"<<endl;
	cout<<"Observed Upper Limit on the ratio R at 95\% CL = "<<rtmp<<endl;
	cout<<"------------------------------------------------------------"<<endl;

	/*
	cms->UseAsimovData();
	double r_asimov = bys.Limit();
	cout<<"Asimov Upper Limit on the ratio R at 95\% CL = "<<r_asimov<<" , "<<endl;
	*/
	if(debug){
		// draw 
		bys.PosteriorPdf();
		TPaveText *pt = SetTPaveText(0.2, 0.7, 0.3, 0.9);
		//start_time=cur_time; cur_time=clock(); cout << "\t TIME_in_BayesianLimit ppdf: "<< (cur_time - start_time)/1000. << " millisec\n";
		DrawEvolution2D pdfr(bys.GetVR(), bys.GetVP(), ";r;Likelihood", "likelihood_pdf", pt);	
		pdfr.setDrawPosteriorPdf(rtmp);
		pdfr.draw();
		pdfr.setDrawPosteriorPdf(rtmp);
		pdfr.save();
	}

	cms->SetSignalScaleFactor(1.);
	if(calcExpectedMeanLimit){
		LimitBands lb(&bys, cms);	
		lb.SetDebug(debug);
		int noutcomes = 1000;
		lb.BysLimitBands(0.05, noutcomes, nexps);

		double rmean, rm1s, rm2s, rp1s, rp2s, rmedian;	
		rmean=lb.GetBysLimitMean();
		rmedian = lb.GetBysLimit(0);
		rm1s=lb.GetBysLimit(-1);rm2s=lb.GetBysLimit(-2);rp1s=lb.GetBysLimit(1);rp2s=lb.GetBysLimit(2);
		cout<<"------------------------------------------------------------"<<endl;
		cout<<"Expected R@95%CL (from -2sigma -1sigma  mean  +1sigma  +2sigma,  median) : "<<endl;
		printf("BANDS  %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n", rm2s, rm1s, rmean, rp1s, rp2s, rmedian);
		cout<<"------------------------------------------------------------"<<endl;


		vector<double> difflimits=lb.GetDifferentialLimitsBys();
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
		pt->AddText("statistical bands");
		string ssave="plot_cump_vs_r95";
		string stitle="; r = #sigma_{95%}/#sigma_{SM}; cumulative probability;";
		PlotXvsCummulativeProb plotRvsP(all_calculated_R95s, cummulativeProbabilities,
				rm1s, rp1s, rm2s, rp2s,ssave, stitle, pt);
		plotRvsP.setGraphDrawOption("CP");
		plotRvsP.draw();

		// show bands without interpolation,  basically using step function
		double GreenBandLow = (1- 0.683)/2.; //1 sigma
		double GreenBandHigh = 1 - (1- 0.683)/2.; //1 sigma
		double YellowBandLow = (1- 0.955)/2.; //2 sigma
		double YellowBandHigh = 1 - (1- 0.955)/2.; //2 sigma
		double rmedian2, rm1s2, rm2s2, rp1s2, rp2s2;	
		bool brmedian2=false, brm1s2=false, brm2s2=false, brp1s2=false, brp2s2=false;	
		for(int i=0; i<cummulativeProbabilities.size(); i++){
			if(cummulativeProbabilities[i]>=GreenBandLow && !brm1s2) { rm1s2 = all_calculated_R95s[i]; brm1s2 = true; } 
			if(cummulativeProbabilities[i]>=GreenBandHigh && !brp1s2) { rp1s2 = all_calculated_R95s[i]; brp1s2 = true; } 
			if(cummulativeProbabilities[i]>=YellowBandLow && !brm2s2) { rm2s2 = all_calculated_R95s[i]; brm2s2 = true; } 
			if(cummulativeProbabilities[i]>=YellowBandHigh && !brp2s2) { rp2s2 = all_calculated_R95s[i]; brp2s2 = true; } 
			if(cummulativeProbabilities[i]>=0.5 && !brmedian2) { rmedian2 = all_calculated_R95s[i]; brmedian2 = true; } 
		}
		cout<<"------------NO INTERPOLATION--------------------------------"<<endl;
		cout<<"BandsNoInterpolation R@95%CL (from -2sigma -1sigma  mean  +1sigma  +2sigma,  median) : "<<endl;
		printf("BANDS2 %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n", rm2s2, rm1s2, rmean, rp1s2, rp2s2, rmedian2);
		cout<<"------------------------------------------------------------"<<endl;

		double rmedian3, rm1s3, rm2s3, rp1s3, rp2s3;	
		GetBandsByLinearInterpolation(all_calculated_R95s, cummulativeProbabilities, rm1s3, rp1s3, rm2s3, rp2s3);
		rmedian3 = GetBandByLinearInterpolation(all_calculated_R95s, cummulativeProbabilities, 0.5);
		cout<<"------------LinearInterpolation--------------------------------"<<endl;
		cout<<"BandsLinearInterpolation R@95%CL (from -2sigma -1sigma  mean  +1sigma  +2sigma,  median) : "<<endl;
		printf("BANDS3 %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n", rm2s3, rm1s3, rmean, rp1s3, rp2s3, rmedian3);
		cout<<"------------------------------------------------------------"<<endl;

		double rmedian4, rm1s4, rm2s4, rp1s4, rp2s4;	
		GetBandsByFeldmanCousins(all_calculated_R95s, cummulativeProbabilities, rm1s4, rp1s4, rm2s4, rp2s4);
		bool brmedian4=false;
		for(int i=0; i<cummulativeProbabilities.size(); i++){
			if(cummulativeProbabilities[i]>=0.5 && !brmedian4) { rmedian4 = all_calculated_R95s[i]; brmedian4 = true; } 
		}
		cout<<"------------FeldmanCousins--------------------------------"<<endl;
		cout<<"BandsByFeldmanCousins R@95%CL (from -2sigma -1sigma  mean  +1sigma  +2sigma,  median) : "<<endl;
		printf("BANDS4 %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n", rm2s4, rm1s4, rmean, rp1s4, rp2s4, rmedian4);
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
									npar++;
									if(argc>=npar+1){
										prior = atoi( argv[npar] );
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
		cout<<"./Bayesian.exe signal err_sig_relative background err_bkg_relative data"<<endl;
		cout<<endl;
		cout<<" or "<<endl;
		cout<<"./Bayesian.exe signal err_sig_relative background err_bkg_relative data calcExpectedMeanLimit=0 pdftypeEs=1 pdftypeEb=1 EsEb_correlated=0 ntoys=100000 seed=1234 debug=0"<<endl;
		cout<<endl<<" detail of the parameters: "<<endl;
		cout<<" signal: 		number of expected signal events " <<endl;
		cout<<" err_sig_relative: 	uncertainty on the signal, in relative fraction " <<endl;
		cout<<" background: 		number of expected background events " <<endl;
		cout<<" err_bkg_relative: 	uncertainty on the background, in relative fraction " <<endl;
		cout<<" data: 			number of observed events"<<endl;
		cout<<" calcExpectedMeanLimit:  calc the mean values of expected limit regardless of the data if 1 " <<endl;
		cout<<" pdftypeEs:		prior pdf on uncertainty of signal,  1 for LogNormal, 2 for TruncatedGaussian"<<endl;
		cout<<" pdftypeEb:		prior pdf on uncertainty of background,  1 for LogNormal, 2 for TruncatedGaussian"<<endl;
		cout<<" EsEb_correlated: 	100\% correlated between Es and Eb if 1,  uncorrelated if 0"<<endl;
		cout<<" ntoys: 			number of toy experiments to avarage out uncertainties "<<endl;
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
	cout<<" seed = "<< seed<<endl;
	cout<<" debuglevel = "<<debug<<endl;
	cout<<"**********************************************"<<endl;
}
