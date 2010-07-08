#include <iostream>
#include "CLsLimit.h"
#include "CRandom.h"
#include "PlotUtilities.h"
#include "FloridaStyle.C"
#include <TError.h>
#include "TSystem.h"
#include "CountingModel.h"
#include "BinnedInterface.h"
#include "BayesianBase.h"
#include "LimitBands.h"
#include "TH1D.h"
using std::cout;
using std::endl;
using namespace lands;
int main(int argc, const char* argv[]){


	// this example will give you instructions on how to run computations for multiple channel combination of shape analysis
	bool calcSignificance=true; 
	bool calcUpperLimits=true;
	int rdmSeed = 1234;
	int debug = 10;
	int ntoys_build_M2lnQ_ForCLsBands = 100000;
	int ntoys_SB_ForSignificance=10000;
	int ntoys_build_M2lnQb_ForSignificance=1000000;
	int ntoys_Bonly_ForLimitBands=1000;
	int ntoys_build_M2lnQ_ForLimitCalc= 10000;
	int ntoys_AverageOutUncertaintiesForBayesian=1000;

	// preparing your input histograms, you can have multiple backgrounds per channel, multipe uncertainties per sample 
	//  in a channel, number of bins and min/max of xaxis should be identical among different samples.
	// you can merge your signal tail as a bin with larger width to reduce the number total bins for speeding up the computation
	TH1D *hsignal[3], *hbackground[3];
	int nbins=100;

	double sig[3]={1.798, 3.481, 2.107};
	double bkg[3] ={ 10.172, 10.342, 12.187};
	double esig[3]={0.110, 0.110, 0.110}; // in relative fraction
	double ebkg[3]={0.28, 0.28, 0.28};
	int pdftypeOfEs[3]={1,1,1}; // 1 for logNormal, 2 for TruncatedGaussian
	int pdftypeOfEb[3]={1,1,1};
	int indexCorrelationEs[3]={1,1,1}; // same number for 100% correlated, 
	int indexCorrelationEb[3]={2,2,2}; // different number for uncorrelated
	CRandom *rdm = new CRandom(rdmSeed);
	for(int c=0; c<3; c++){
		TString s = (TString)c;
		hsignal[c]     = new TH1D("signal"+c,"the expected signal"+c,nbins,-5,5);
		hbackground[c] = new TH1D("background1"+c,"The expected background"+c,nbins,-5,5);
		// Fill in templates, i.e. expected signal and background shape, normalized to their expected number of events at certain intL
		for(int n=0; n<100000; n++){
			hsignal[c]->Fill(rdm->Gaus(0,1), sig[c]/100000.);	
			hbackground[c]->Fill(rdm->Rndm() * 10 -5, bkg[c]/100000.);	
		}
	}


	//  constructing model
	CountingModel *cms=new CountingModel();
	cms->SetRdm(rdm);
	BinnedInterface *bi = new BinnedInterface(cms);
	for(int c=0; c<3; c++) {
		bi->AddChannel(hsignal[c], hbackground[c]);
		// bi->AddUncertainty(int index_channel, int index_sample, double error_in_relative_fraction, int pdf_type, int index_correlation);
		// index_sample: in a channel,  0 means signal, 1-5 means background 1-5, currently we set maximum number of categories of backgrounds at 5 
		// error_in_relative_fraction: assume that for this paticular uncertainty source, all bins in this sample has same relative error scale
		//				if not, then,  you may use bi->AddUncertainty(int index_channel, int index_sample, TH1* h_error_in_relative_fraction, int pdf_type, int index_correlation);
		//				That is good for uncertainty from different shapes
		// pdf_type:  1 for LogNormal, 2 for TruncatedGaussian
		// index_correlation:   among channels and among samples,  same number for 100% correlation, different number for 0% correlation
		bi->AddUncertainty(c, 0, esig[c], pdftypeOfEs[c], indexCorrelationEs[c]);
		bi->AddUncertainty(c, 1, ebkg[c], pdftypeOfEb[c], indexCorrelationEb[c]);
	}

	// We don't consider systematics errors by default, you can switch it as you like
	cms->SetUseSystematicErrors(true);


	//  run and test ....
	CLsBase frequentist;
	frequentist.SetDebug(debug);

	CLsLimit clsr95;
	clsr95.SetDebug(debug);

	// to store results of the statistical bands 
	double rmean, rm1s,rm2s,rp1s,rp2s;  //mean, -1sigma, -2sigma, +1sigma, +2sigma


	clsr95.DoingStatisticalBandsForCLs(cms, &frequentist, ntoys_build_M2lnQ_ForCLsBands); // ntoys_to_build_-2lnQ
	rmean=clsr95.CLs_mean(); rm2s=clsr95.CLs_sigma(-2); rm1s=clsr95.CLs_sigma(-1); 	rp1s=clsr95.CLs_sigma(1); rp2s=clsr95.CLs_sigma(2);
	cout<<"------------------------------------------------------------"<<endl;
	cout<<"Expected CLs (from -2sigma -1sigma  mean  +1sigma  +2sigma) : "<<endl;
	printf("----------- %10.5f %10.5f %10.5f %10.5f %10.5f\n", rm2s, rm1s, rmean, rp1s, rp2s); 
	cout<<"------------------------------------------------------------"<<endl;

	vector<double> vm2logQ_b, vm2logQ_b_prob;
	SortAndCumulative(frequentist.Get_m2logQ_b(), vm2logQ_b, vm2logQ_b_prob, 1);// sort it by  decreased order // this frequentist is run done by previous step
	GetBandsByLinearInterpolation(vm2logQ_b,vm2logQ_b_prob, rm1s, rp1s, rm2s, rp2s );
	rmean=GetMeanOfSortedXwithProb(vm2logQ_b, vm2logQ_b_prob);
	cout<<"------------------------------------------------------------"<<endl;
	cout<<"Expected -2lnQ for background only (from -2sigma -1sigma  mean  +1sigma  +2sigma) : "<<endl;
	printf("----------- %10.5f %10.5f %10.5f %10.5f %10.5f\n", rm2s, rm1s, rmean, rp1s, rp2s); 
	cout<<"------------------------------------------------------------"<<endl;

	vector<double> vm2logQ_sb, vm2logQ_sb_prob;
	SortAndCumulative(frequentist.Get_m2logQ_sb(), vm2logQ_sb, vm2logQ_sb_prob, 1); // sort it by  decreased order
	GetBandsByLinearInterpolation(vm2logQ_sb,vm2logQ_sb_prob, rm1s, rp1s, rm2s, rp2s );
	rmean=GetMeanOfSortedXwithProb(vm2logQ_sb, vm2logQ_sb_prob);
	cout<<"------------------------------------------------------------"<<endl;
	cout<<"Expected -2lnQ for bkg+signal (from -2sigma -1sigma  mean  +1sigma  +2sigma) : "<<endl;
	printf("----------- %10.5f %10.5f %10.5f %10.5f %10.5f\n", rm2s, rm1s, rmean, rp1s, rp2s); 
	cout<<"------------------------------------------------------------"<<endl;


	if(calcSignificance){
		cout<<"\t Start calculation of significance ... "<<endl;
		vector<double> vsignificance, vsignificance_cp;
		rmean = frequentist.SignificanceComputation(ntoys_SB_ForSignificance, ntoys_build_M2lnQb_ForSignificance, vsignificance, vsignificance_cp );//mean value
		GetBandsByLinearInterpolation(vsignificance,vsignificance_cp, rm1s, rp1s, rm2s, rp2s );//68% and 95% bands
		double rmedian=GetBandByLinearInterpolation(vsignificance, vsignificance_cp, 0.5);//median value
		cout<<"------------------------------------------------------------"<<endl;
		cout<<"Expected significance (from -2sigma -1sigma  median  +1sigma  +2sigma  mean) : "<<endl;
		printf("----------- %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n", rm2s, rm1s, rmedian, rp1s, rp2s, rmean); 
		cout<<"------------------------------------------------------------"<<endl;
	}

	BayesianBase bys;
	bys.SetDebug(debug);
	if(calcUpperLimits){
		LimitBands lb(&clsr95, &frequentist, &bys, cms);	
		//Bands(1-C.L., n_sets_outcomes, doCLsLimitBands, ntoys_buildM2lnQ, doBysLimitBands, ntoys_averageOutUncertainties) 
		lb.Bands(0.05, ntoys_Bonly_ForLimitBands, true, ntoys_build_M2lnQ_ForLimitCalc, true, ntoys_AverageOutUncertaintiesForBayesian);

		rmean=lb.GetCLsLimitMean();
		rm1s=lb.GetCLsLimit(-1);rm2s=lb.GetCLsLimit(-2);rp1s=lb.GetCLsLimit(1);rp2s=lb.GetCLsLimit(2);
		cout<<"------------------------------------------------------------"<<endl;
		cout<<"Expected CL95% UpperLimit (Modified Frequentist approach) (from -2sigma -1sigma  mean  +1sigma  +2sigma) : "<<endl;
		printf("----------- %10.5f %10.5f %10.5f %10.5f %10.5f\n", rm2s, rm1s, rmean, rp1s, rp2s); 
		cout<<"------------------------------------------------------------"<<endl;

		rmean=lb.GetBysLimitMean();
		rm1s=lb.GetBysLimit(-1);rm2s=lb.GetBysLimit(-2);rp1s=lb.GetBysLimit(1);rp2s=lb.GetBysLimit(2);
		cout<<"------------------------------------------------------------"<<endl;
		cout<<"Expected CL95% UpperLimit (Bayesian approach) (from -2sigma -1sigma  mean  +1sigma  +2sigma) : "<<endl;
		printf("----------- %10.5f %10.5f %10.5f %10.5f %10.5f\n", rm2s, rm1s, rmean, rp1s, rp2s); 
		cout<<"------------------------------------------------------------"<<endl;
	}
	return 1;
}
