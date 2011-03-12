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
#include "TGraph.h"
using std::cout;
using std::endl;
using namespace lands;
void processParameters(int argc, const char* argv[]);
int debug=0; int nexps=100000; double s=0;  double b=0; double s_err = 0; double b_err = 0; int d=0;
int seed =1234; int pdftypeEs = 1; int pdftypeEb = 1; int EsEb_correlated = 0; int calcExpectedMeanLimit=0; 
int testStatistics = 1, rule = 1; // default is CLs
int asimov = -1; // -1 for using the input dataset whatever you provide,  0 for using asimov bkg only dataset, 1 for asimov sig+bkg

// current for FeldmanCousins
bool bAdaptiveSampling = true;
bool bQuickEstimateInitialLimit = true;
double toysFactor_ = 1.;
double rAbsAccuracy_ = 0.01, rRelAccuracy_ = 0.01;
//double rAbsAccuracy_ = 0.1;
//double rRelAccuracy_ = 0.05;
bool lowerLimit_ = false;
double intialRmin = 1., intialRmax = 21;// only when bQuickEstimateInitialLimit==false

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

	//cms->SetAllowNegativeSignalStrength(false);

	// initialize the calculator
	CLsBase frequentist;
	frequentist.SetDebug(debug);
	frequentist.SetRdm(rdm);
	frequentist.SetTestStatistics(testStatistics);
	cms_global= cms;
	//vdata_global=cms->Get_v_data();

	double tmp;
	vdata_global=cms->Get_v_data();


	if(rule==1 || rule==2)frequentist.BuildM2lnQ(cms,nexps);
	else frequentist.BuildM2lnQ(cms, 100);
	double errs, errb, errsb;
	double cls = frequentist.CLs(errs);
	double clsb = frequentist.CLsb(errsb);
	double clb = frequentist.CLb(errb);
	cout<<"------------------------------------------------------------"<<endl;
	cout<<"Observed CLs = "<<cls<<" +/- "<<errs<<endl;
	cout<<"Observed CLsb = "<<clsb<<" +/- "<<errsb<<endl;
	cout<<"Observed CLb = "<<clb<<" +/- "<<errb<<endl;
	cout<<"------------------------------------------------------------"<<endl;

	DrawPdfM2logQ pdfM2logQ(frequentist.Get_m2logQ_sb(),frequentist.Get_m2logQ_b(), frequentist.Get_m2lnQ_data(), 
			"-2lnQ on data", "; -2lnQ; entries", "lnq", pt);
	pdfM2logQ.draw();
	double m2lnqdata = frequentist.Get_m2lnQ_data();
	//cout<<m2lnqdata<<endl;
	vector<double> vsb = frequentist.Get_m2logQ_sb();
	for(int i=0; i<vsb.size(); i++){
	//	if(vsb[i]>=m2lnqdata)cout<<i<<" "<<vsb[i]<<endl;
	}
	//cout<<m2lnqdata<<endl;

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

	if(rule==3){

		// FeldmanCousins must use specified test statisics:  profile likelihood ratio, which allow mu hat > probed mu
		frequentist.SetTestStatistics(31);

		cms->SetAllowNegativeSignalStrength(false);
		clsr95.SetAdaptiveSampling(bAdaptiveSampling);
		bool lowerLimit_ = false;
	
		double r95_fc;

		double fcMid = 0, fcErr = 0; 
		double rmin = intialRmin,  rmax = intialRmax;

		if(bQuickEstimateInitialLimit){
			BayesianBase bys(cms, 0.05, 1.e-2);
			bys.SetNumToys(1000);
			rtmp= bys.Limit();
			rmin = rtmp/10.,  rmax = rtmp*2;
		}

		int nsteps = 10;
		vector< vector<double> > vv_all;
		vector< vector<double> > vv;
		// the following iterating algothm comes from http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/HiggsAnalysis/CombinedLimit/src/FeldmanCousins.cc
		do { 
			if (debug) std::cout << "scan in range [" << rmin << ", " << rmax << "]" << std::endl;
			if(fcMid!=0){
				double bin = (rmax-rmin)/10.;
				rmin = rmin+0.5*bin;
				rmax = rmax-0.5*bin;
			}
			r95_fc = clsr95.FeldmanCousins(cms, rmin, rmax, &frequentist, nexps, nsteps);
			vv = clsr95.GetFCconstruction();
			for(int i=0; i<vv.size(); i++){
				vv_all.push_back(vv[i]);
			}
			int found = -1; 
			if (lowerLimit_) {
				for(int i=0; i<nsteps; ++i){
					int n = vv[i].size();
					if(vv[i][n-1]>vv[i][n-2]) found=i;
					else break;
				}
			} else {
				for(int i=0; i<nsteps; ++i){
					int n = vv[i].size();
					bool inside = (vv[i][n-1]<=vv[i][n-2]) ; // important ...
					if (inside) found = i;
					else if (found != -1) break;
				}
				if (found == -1) {
					std::cout << "Points are either all inside or all outside the bound." << std::endl;
					return false;
				}
			}
			double fcBefore = (found > -1 ? (rmin+(rmax-rmin)/double(nsteps)*found) : rmin);
			double fcAfter  = (found < nsteps-1 ? (rmin+(rmax-rmin)/double(nsteps)*(found+1)) : rmax);
			fcMid = 0.5*(fcAfter+fcBefore);
			fcErr = 0.5*(fcAfter-fcBefore);
			if (debug) std::cout << "  would be r < " << fcMid << " +/- "<<fcErr << std::endl;
			rmin=(std::max(rmin, fcMid-3*fcErr)); 
			rmax=(std::min(rmax, fcMid+3*fcErr));
			if (fcErr < 4*std::max(rAbsAccuracy_, rRelAccuracy_ * fcMid)) { // make last scan more precise
				clsr95.SetAdditionalNToysFactor(4*toysFactor_);
			}
		} while (fcErr > std::max(rAbsAccuracy_, rRelAccuracy_ * fcMid));

		r95_fc = fcMid;
		if (1) {
			std::cout << "\n -- FeldmanCousins++ -- \n";
			std::cout << "Limit: r " << (lowerLimit_ ? "> " : "< ") << fcMid << " +/- "<<fcErr << "\n";
		}



		// for plots...
		TCanvas *can = new TCanvas("c","c");	
		TString ssave = "fc";
		TString stmp = "hframe"; stmp += ssave;
		double xmax = 0;
		for(int i=0; i<vv_all.size(); i++){
			int n = vv_all[i].size();
			if(vv_all[i][n-2]>xmax) xmax= vv_all[i][n-2];
		}
		TH1F *hframe= new TH1F(stmp, "; q_{#mu}; #mu", 1000, 0, 2*xmax);
		hframe->SetMinimum(rtmp/10.);
		hframe->SetMaximum(rtmp*2.);
		hframe->SetStats(0);
		hframe->SetFillStyle(1);
		hframe->Draw(" ");

		float *r_tested = new float[vv_all.size()];
		float *q_data = new float[vv_all.size()];

		float q95=0;
		for(int i=0; i<vv_all.size(); i++){
			int n = vv_all[i].size();
			double x[2]={0, vv_all[i][n-2]};	// q_up	
			double y[2]={vv_all[i][n-3], vv_all[i][n-3]}; // r
			TGraph *gr = new TGraph(2, x, y);
			gr->Draw("l same");

			r_tested[i]=vv_all[i][n-3];
			q_data[i]=vv_all[i][n-1];

			//if(r95_fc==vv_all[i][n-3]) q95 = q_data[i];
			if(fabs(r95_fc-r_tested[i])/r95_fc<0.00001) q95 = q_data[i];
		}

		TGraph gr(vv_all.size(), q_data, r_tested);
		gr.Sort(&TGraph::CompareY);
		gr.Draw("l same");

		TArrow *arrow95 = new TArrow(0, r95_fc, q95, r95_fc, 0.03, "|>");
		arrow95->SetLineWidth(3);
		arrow95->SetLineColor(kRed);
		//arrow95->SetFillColor(kRed);
		arrow95->SetFillStyle(0);
		arrow95->Draw();
		Save(can, "fc");

		cout <<"95\% CL upper limit by FC: "<<r95_fc<<",   use bys: "<<rtmp<<endl;

		delete r_tested;
		delete q_data;
		return 1;
	}



	rtmp = clsr95.LimitOnSignalScaleFactor(cms, &frequentist,nexps);

	if(debug){
		DrawEvolution2D d2d(clsr95.GetvTestedScaleFactors(), clsr95.GetvTestedCLs(), "; r ; CLs", "r_vs_cls", pt);
		d2d.draw();
		d2d.save();

	}

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
		double rmedian = lb.GetBysLimit(0);
		rm1s=lb.GetCLsLimit(-1);rm2s=lb.GetCLsLimit(-2);rp1s=lb.GetCLsLimit(1);rp2s=lb.GetCLsLimit(2);
		cout<<"------------------------------------------------------------"<<endl;
		cout<<"Expected R@95%CL (from -2sigma -1sigma  mean  +1sigma  +2sigma) : "<<endl;
		printf("--------- %10.5f %10.5f %10.5f %10.5f %10.5f\n", rm2s, rm1s, rmean, rp1s, rp2s);
		printf("BANDS  %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n", rm2s, rm1s, rmean, rp1s, rp2s, rmedian);
		cout<<"------------------------------------------------------------"<<endl;
	
		TString ts(fileName); ts+="_clslimits";
		FillTree(ts, lb.GetDifferentialLimitsCLs());

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
		cout<<"./CLs_dataCard.exe inputFileName"<<endl;
		cout<<endl;
		cout<<" or "<<endl;
		cout<<"./CLs_dataCard.exe inputFileName calcExpectedMeanLimit=0 ntoys=100000 seed=1234 debug=0 testStatistics=1 rule=1"<<endl;
		cout<<endl<<" detail of the parameters: "<<endl;
		cout<<" inputFileName:          The input data card designed by Andrey Korytov "<<endl;
		cout<<" calcExpectedMeanLimit:  calc the mean values of expected limit regardless of the data if 1 " <<endl;
		cout<<" ntoys: 			number of toy experiments to build -2lnQ "<<endl;
		cout<<" seed:			seed for random number generator "<<endl;
		cout<<" debug: 			debug level "<<endl;
		cout<<" testStatistics:		1 for Q_LEP, 2 for Q_TEV, 4 for only profile mu"<<endl;
		cout<<"                         3 for Q_ATLAS, 31 for Q_ATLAS but allowing mu_hat>mu"<<endl;
		cout<<" rule:                   1 for CLs,  2 for CLsb, 3 for FeldmanCousins"<<endl;
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
	else if(testStatistics==31) cout<<" ATLAS type but allowing mu_hat>mu";
	else if(testStatistics==4) cout<<" only profile mu";
	else cout<<" LEP type";
	cout<<endl;

	cout<<" rule = "<< rule<<", is ";
	if(rule==2) cout<<" CLsb";
	else if(rule==1) cout<<" CLs";
	else if(rule==3) cout<<" FeldmanCousins";
	cout<<endl;

	if(asimov==0) cout<<" Using Asimov bkg only dataset"<<endl;
	if(asimov==1) cout<<" Using Asimov signal + bkg  dataset"<<endl;

	cout<<"**********************************************"<<endl;
}
