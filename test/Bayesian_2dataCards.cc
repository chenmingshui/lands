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
const char * fileName, *fileName2;
int main(int argc, const char* argv[]){
	processParameters(argc, argv);

//	FloridaStyle();

	CRandom *rdm = new CRandom(seed);  //initilize a random generator

	// set up couting model
	CountingModel *cms=new CountingModel();
	CountingModel *cms1=new CountingModel();
	CountingModel *cms2=new CountingModel();

	ConfigureModel(cms1, fileName); 
	ConfigureModel(cms2, fileName2); 
	cms1->SetUseSystematicErrors(true);
	cms2->SetUseSystematicErrors(true);
	cout<<"before combine"<<endl;
	CountingModel cms3 = CombineModels(cms1, cms2);
	cout<<"after combine"<<endl;
	cms = &cms3;

	// everything else for cms need to be reset, because now cms pointing to cms3 

	cms->SetRdm(rdm);
	cout<<"&cms3"<<endl;
	cms->SetUseSystematicErrors(true);
	cout<<"good"<<endl;

	cms->UseAsimovData(asimov);

	cms->Print(debug);
	cms->RemoveChannelsWithExpectedSignal0orBkg0();

	cms->Print(debug);
	
	//cms->SetDebug(10);
	BayesianBase bys(cms, 0.05, 1.e-3);
	bys.SetNumToys(nexps);
	bys.SetDebug(debug);
	if(XSprior==20)bys.SetCrossSectionPrior(corr);
	if(XSprior==10)bys.SetCrossSectionPrior(flat);
	if(XSprior==30)bys.SetCrossSectionPrior(prior_1overSqrtS);
	double rtmp;
	rtmp = bys.Limit();
	cout<<"------------------------------------------------------------"<<endl;
	cout<<"Observed Upper Limit on the ratio R at 95\% CL = "<<rtmp<<endl;
	cout<<"------------------------------------------------------------"<<endl;

	// draw 

	if(debug)	{
		bys.PosteriorPdf();
		TPaveText *pt = SetTPaveText(0.2, 0.7, 0.3, 0.9);
		//start_time=cur_time; cur_time=clock(); cout << "\t TIME_in_BayesianLimit ppdf: "<< (cur_time - start_time)/1000. << " millisec\n";
		DrawEvolution2D pdfr(bys.GetVR(), bys.GetVP(), ";r;Likelihood", "likelihood_pdf", pt);	
		pdfr.setDrawPosteriorPdf(rtmp);
		pdfr.draw();
		pdfr.getGraph()->SetMarkerSize(0.01);
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
		rm1s=lb.GetBysLimit(-1);rm2s=lb.GetBysLimit(-2);rp1s=lb.GetBysLimit(1);rp2s=lb.GetBysLimit(2);
		rmedian = lb.GetBysLimit(0);
		// plot the distribution of limits of noutcomes
		// and the cummulative pdf
		TString ts(fileName); ts+="_and_"; ts+=fileName2; ts+="_byslimits";
		FillTree(ts, lb.GetDifferentialLimitsBys());

		// default is with interpolation using Fermi function
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
		//pt->AddText("statistical bands");
		string ssave="plot_cump_vs_r95";
		string stitle="; r = #sigma_{95%}/#sigma_{SM}; cumulative probability;";
		PlotXvsCummulativeProb plotRvsP(all_calculated_R95s, cummulativeProbabilities,
				rm1s, rp1s, rm2s, rp2s,ssave, stitle, pt);
		plotRvsP.draw();
		plotRvsP.getGraph()->SetMarkerSize(1);
		plotRvsP.save();
		
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

		ssave="plot_cump_vs_r95_step";
		stitle="; r = #sigma_{95%}/#sigma_{SM}; cumulative probability;";
		PlotXvsCummulativeProb plotRvsP2(all_calculated_R95s, cummulativeProbabilities,
				rm1s2, rp1s2, rm2s2, rp2s2,ssave, stitle, pt);
		plotRvsP2.setGraphDrawOption("P");
		plotRvsP2.draw();
		plotRvsP2.getCanvas()->cd();
		plotRvsP2.getGraph()->SetMarkerSize(1);
		int ntmp = all_calculated_R95s.size();
		double *xtmp = new double[ntmp];
		double *ytmp = new double[ntmp];
		for(int i=0; i<ntmp; i++){
			xtmp[i]=all_calculated_R95s[i];
			ytmp[i]=cummulativeProbabilities[i];
		}
		plotRvsP2.getGraph()->PaintGrapHist(ntmp, xtmp, ytmp, "hfsame");
		plotRvsP2.save();

		double rmedian3, rm1s3, rm2s3, rp1s3, rp2s3;	
		GetBandsByLinearInterpolation(all_calculated_R95s, cummulativeProbabilities, rm1s3, rp1s3, rm2s3, rp2s3);
		rmedian3 = GetBandByLinearInterpolation(all_calculated_R95s, cummulativeProbabilities, 0.5);
		cout<<"------------LinearInterpolation--------------------------------"<<endl;
		cout<<"BandsLinearInterpolation R@95%CL (from -2sigma -1sigma  mean  +1sigma  +2sigma,  median) : "<<endl;
		printf("BANDS3 %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n", rm2s3, rm1s3, rmean, rp1s3, rp2s3, rmedian3);
		cout<<"------------------------------------------------------------"<<endl;

		ssave="plot_cump_vs_r95_linear";
		stitle="; r = #sigma_{95%}/#sigma_{SM}; cumulative probability;";
		PlotXvsCummulativeProb plotRvsP3(all_calculated_R95s, cummulativeProbabilities,
				rm1s3, rp1s3, rm2s3, rp2s3,ssave, stitle, pt);
		plotRvsP3.setGraphDrawOption("LP");
		plotRvsP3.draw();
		plotRvsP3.getGraph()->SetMarkerSize(1);
		plotRvsP3.save();

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

		ssave="plot_cump_vs_r95_FC";
		stitle="; r = #sigma_{95%}/#sigma_{SM}; cumulative probability;";
		PlotXvsCummulativeProb plotRvsP4(all_calculated_R95s, cummulativeProbabilities,
				rm1s4, rp1s4, rm2s4, rp2s4,ssave, stitle, pt);
		plotRvsP4.setGraphDrawOption("LP");
		plotRvsP4.draw();
		plotRvsP4.getGraph()->SetMarkerSize(1);
		plotRvsP4.save();
	}

	return 1;
}
void processParameters(int argc, const char* argv[]){
	int npar=1;
	if(argc>=3){
		fileName = argv[1];
		npar++;
		fileName2 = argv[2];
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
							if(XSprior !=10 and XSprior!=20 and XSprior!=30) {
								cout<<"prior on the cross section should be \"10\" or \"20\" or \"30\""<<endl;
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
		cout<<"./Bayesian_dataCard.exe inputFileName1 inputFileName2"<<endl;
		cout<<endl;
		cout<<" or "<<endl;
		cout<<"./Bayesian_dataCard.exe inputFileName1 inputFileName2 calcExpectedMeanLimit=0 ntoys=100000 seed=1234 debug=0 XSprior=10"<<endl;
		cout<<endl<<" detail of the parameters: "<<endl;
		cout<<" inputFileName1/2:          The input data card designed by Andrey Korytov "<<endl;
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
	cout<<" Input data cards are:   "<<fileName<<" and "<<fileName2<<endl; 
	cout<<" prior on the cross section ratio = ";
		if(XSprior==10) cout<<"flat"<<endl;
		if(XSprior==20) cout<<"correlated"<<endl;
		if(XSprior==30) cout<<"1/sqrt(r)"<<endl;
	cout<<(calcExpectedMeanLimit==0?" NOT":"")<<" do calculation of mean limit"<<endl;
	cout<<" seed = "<< seed<<endl;
	cout<<" debuglevel = "<<debug<<endl;

	if(asimov==0) cout<<" Using Asimov bkg only dataset"<<endl;
	if(asimov==1) cout<<" Using Asimov signal + bkg  dataset"<<endl;

	cout<<"**********************************************"<<endl;
}
