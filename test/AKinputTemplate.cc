#include <iostream>
#include <vector>
#include "CLsLimit.h"
#include "CRandom.h"
#include "PlotUtilities.h"
#include "CountingModel.h"
#include "BayesianBase.h"
#include "LimitBands.h"
#include "UtilsROOT.h"
#include "TMath.h"

#include "TGraph.h"

#include <time.h> // upto micro second

#include "TString.h"
#include "TObjString.h"
#include "TObject.h"
#include "TObjArray.h"
#include "TIterator.h"
#include <fstream>

using std::cout;
using std::endl;
using namespace lands;

double StandarDeviation(int n, const double* d){
	double mean = TMath::Mean(n, d);
	double v = 0;
	for(int i=0; i<n; i++){
		v+=( (d[i]-mean)*(d[i]-mean) );
	}
	return sqrt(v/(n));
}
int main(int argc, const char* argv[]){

	const char* fileName =  "/data/Projects/Statistics/LandS/test/ATLAS_CMS/cmsinput/2010.06.22_hww_mh140A_3ch_0jets.txt";
	if(argc>=2) fileName = argv[1];


	clock_t start_time=clock(), cur_time=clock();


	//  configure here what do you want to do:   how many toys to generate, etc....
	bool calcSignificance=false; 
	bool calcUpperLimits=true;
	int rdmSeed = 1234;
	int debug = 0;
	int ntoys_build_M2lnQ_ForCLsBands = 100000;
	int ntoys_SB_ForSignificance=10000;
	int ntoys_build_M2lnQb_ForSignificance=100000000;
	int ntoys_Bonly_ForLimitBands=1000;
	int ntoys_build_M2lnQ_ForLimitCalc= 100000;
	int ntoys_AverageOutUncertaintiesForBayesian=20000;

	CRandom *rdm = new CRandom(rdmSeed);
	CountingModel* cms=new CountingModel();
	cms->SetRdm(rdm);

	// configure the model from a template input file
	ConfigureModel(cms, fileName); 


	// set to be using nuisance parameters
	cms->SetUseSystematicErrors(true);
	cms->Print();

	TPaveText *pt = SetTPaveText(0.2, 0.7, 0.3, 0.9);
	double rtmp;

	//  run and test ....
	CLsBase frequentist;
	frequentist.SetDebug(debug);

	frequentist.SetModel(cms);


	// precisely calculate the significance without systematics
	// ----   quick test to get significance without uncertainties and for small number of channels
	/* 
	int nmax = 100;
	double totalbkg[3]={0,0,0};
	for(int i=0; i<6; i++) {
		for(int j=0; j<3; j++){
			totalbkg[j]+=bkg[j][i];
		}
	}
	   frequentist.BuildM2lnQ(cms,1);
	   double lnq_data = frequentist.Get_m2lnQ_data();

	   TCanvas *ctemp = new TCanvas("ctemp", "ctemp");	
	   TH1F *htemp = new TH1F("htemp",";-2lnQ; ", 12000, -60, 60);

	   double prob_data = 0;
	   double minus_prob_data = 0;
	   double lnqtmp = 0, p=0;
	   for(int i=0; i<100; i++){
	   cout<<"CountLoop "<< i<<endl;
	   for(int j=0; j<100; j++)
	   for(int k=0; k<100; k++){
	   
	   //reset data numbers 
	   cms->AddObservedData(0, i);
	   cms->AddObservedData(1, j);
	   cms->AddObservedData(2, k);
	   frequentist.BuildM2lnQ(cms,1);
	   lnqtmp = frequentist.Get_m2lnQ_data() ;
	   p = TMath::Poisson(i, totalbkg[0]) * TMath::Poisson(j, totalbkg[1]) * TMath::Poisson(k, totalbkg[2]) ; 
	   htemp  -> Fill (lnqtmp, p);
	   if( lnq_data > lnqtmp )
	   prob_data += p ;
	   else 
	   minus_prob_data += p;
	   }
	   }

	   ctemp->SetLogy(1);
	   htemp->Draw();
	   Save(ctemp, "analytical_lnq");

	   cout<<"p_value = "<<prob_data <<",  1-p_value "<<minus_prob_data<<endl;
	   cout<<"significance = "<<Significance(prob_data)<<endl;


	   return 1;
	   */



	CLsLimit clsr95;
	clsr95.SetDebug(debug);

	int nexps = ntoys_build_M2lnQ_ForCLsBands;
	frequentist.BuildM2lnQ(cms,nexps);
	double cls = frequentist.CLs();
	double clb = frequentist.CLb();
	double clsb = frequentist.CLsb();
	cout<<"------------------------------------------------------------"<<endl;
	cout<<"Observed CLsb = "<<clsb<<endl;
	cout<<"Observed CLb  = "<<clb<<endl;
	cout<<"Observed CLs  = "<<cls<<endl;
	cout<<"------------------------------------------------------------"<<endl;

	//return 1;
	DrawPdfM2logQ pdfM2logQ(frequentist.Get_m2logQ_sb(),frequentist.Get_m2logQ_b(), frequentist.Get_m2lnQ_data(), 
			"-2lnQ on data", "; -2lnQ; entries", "lnq", pt);
	pdfM2logQ.draw();

	start_time=cur_time; cur_time=clock(); cout << "\t\t\t  " << (cur_time - start_time)/1000. << " millisec\n"; 
	clsr95.SetAlpha(0.05);
	rtmp = clsr95.LimitOnSignalScaleFactor(cms, &frequentist,nexps);
	//rtmp = clsr95.LimitOnSignalScaleFactor(cms, 0.2, 0.4, &frequentist,nexps );
	cout<<"------------------------------------------------------------"<<endl;
	cout<<"Observed Upper Limit on the ratio R at 95\% CL = "<<rtmp<<endl;
	cout<<"------------------------------------------------------------"<<endl;

	start_time=cur_time; cur_time=clock(); cout << "\t\t\t CLs limit " << (cur_time - start_time)/1000. << " millisec\n"; 

	vector<double> vr = clsr95.GetvTestedScaleFactors();	
	DrawEvolution2D drawRvsCL(vr, clsr95.GetvTestedCLs(), "; r = #sigma/#sigma_{SM};CLs;", "rvscls", pt);
	drawRvsCL.draw();


	start_time=cur_time; cur_time=clock(); cout << "\t\t\t " << (cur_time - start_time)/1000. << " millisec\n"; 

	// bayesian, need reset the observed data
	/*   
	// reset signal, bkg numbers
	cms->SetSignalScaleFactor(1);
	cms->Print();
	BayesianBase bys(cms, 0.05, 1.e-3);
	bys.SetNumToys(20000);
	bys.SetDebug(debug);
	//bys.SetCrossSectionPrior(corr);
	bys.SetCrossSectionPrior(flat);
	rtmp = bys.Limit();
	cout<<"------------------------------------------------------------"<<endl;
	cout<<"Observed Upper Limit on the ratio R at 95\% CL = "<<rtmp<<endl;
	cout<<"------------------------------------------------------------"<<endl;

	start_time=cur_time; cur_time=clock(); cout << "\t\t\t Bayesian limit " << (cur_time - start_time)/1000. << " millisec\n"; 
	double mean =0 ;
	cout<<" 20 trials on Bayesian limits and its standard deviation is =  "<<bys.ErrorOnR_DueToFiniteToys(20, mean)<<endl;
	cout<<" its mean = "<<mean<<endl;

	start_time=cur_time; cur_time=clock(); cout << "\t\t\t MEAN Bayesian limit " << (cur_time - start_time)/1000. << " millisec\n"; 
	// draw 

	double err_on_r = bys.PosteriorPdf();
	cout<<"err_on_bayesian_r (due to Alpha precision) = "<< err_on_r <<endl;
	//start_time=cur_time; cur_time=clock(); cout << "\t TIME_in_BayesianLimit ppdf: "<< (cur_time - start_time)/1000. << " millisec\n";
	DrawEvolution2D pdfr(bys.GetVR(), bys.GetVP(), ";r;Likelihood", "likelihood_pdf", pt);	
	pdfr.draw();
	pdfr.getGraph()->Draw("AL");
	pdfr.getLine()->Delete();
	pdfr.save();
	*/


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

	if(calcUpperLimits && 0){
		BayesianBase bys;
		bys.SetDebug(debug);
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
