#include <iostream>
#include <vector>
#include "CLsLimit.h"
#include "CRandom.h"
#include "PlotUtilities.h"
#include "CountingModel.h"
#include "BayesianBase.h"
#include "LimitBands.h"
#include "UtilsROOT.h"

using std::cout;
using std::endl;
using namespace lands;

int main(int argc, const char* argv[]){


	// this example will give you instructions on how to run computations for multiple channel combination
	bool calcSignificance=true; 
	bool calcUpperLimits=true;
	int rdmSeed = 1234;
	int debug = 0;
	int ntoys_build_M2lnQ_ForCLsBands = 100000;
	int ntoys_SB_ForSignificance=10000;
	int ntoys_build_M2lnQb_ForSignificance=10000000;
	int ntoys_Bonly_ForLimitBands=1000;
	int ntoys_build_M2lnQ_ForLimitCalc= 100000;
	int ntoys_AverageOutUncertaintiesForBayesian=20000;

	// inputs
	double sig[3]={1.798, 3.481, 2.107};
	double bkg[3] ={ 10.172, 10.342, 12.187};
	double esig[3]={0.110, 0.110, 0.110}; // in relative fraction
	double ebkg[3]={0.28, 0.28, 0.28};
	int pdftypeOfEs[3]={1,1,1}; // 1 for logNormal, 2 for TruncatedGaussian
	int pdftypeOfEb[3]={1,1,1};
	int indexCorrelationEs[3]={1,1,1}; // same number for 100% correlated, 
	int indexCorrelationEb[3]={2,2,2}; // different number for uncorrelated

	// constructing model
	CountingModel* cms=new CountingModel();
	CRandom *rdm = new CRandom(rdmSeed);
	cms->SetRdm(rdm);
	for(int c=0; c<3; c++){
		cms->AddChannel(sig[c], bkg[c]);
		cms->AddUncertainty(c, 0, esig[c], pdftypeOfEs[c], indexCorrelationEs[c]);
		cms->AddUncertainty(c, 1, ebkg[c], pdftypeOfEb[c], indexCorrelationEb[c]);
	}	
	cms->SetUseSystematicErrors(true);
	if(debug)cms->Print();




	
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
