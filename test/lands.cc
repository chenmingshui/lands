#include <TString.h>
#include <map>
#include <stdio.h>
#include <vector>
#include <iostream>
#include <TSystem.h>
#include "CountingModel.h"
#include "BayesianBase.h"
#include "LimitBands.h"
#include "CRandom.h"
#include "UtilsROOT.h"
#include "PlotUtilities.h"
#include "CLsLimit.h"
#include "SignificanceBands.h"
//#include <cmath>

using namespace std;
using namespace lands;

typedef std::map<TString, vector<TString> > TMap;
typedef std::pair<TString, vector<TString> > TPair;
TMap options;
void processParameters(int argc, const char* argv[]);
void PrintHelpMessage();

// parameters :
vector<TString> datacards; // enable multiple cards combination, and allow "* ?" in the input
TString method; // ProfiledLikelihood, Bayesian, Hybrid, FeldmanCousins 
TString jobname; // default constructed as "fist_card_name"+"method"
TString dataset; // default: data_obs;  asimov_b, asimov_sb
int debug; // default 0
bool doExpectation; // default false; 
int toys; // number of toys to do expectation, e.g. bands and mean/median limits.   default 1000,  you must also turn on doExpectation
int toysHybrid; // number of toys used to build CLsb, CLb, CLs, as well as number of toys per point in FeldmanCousins construction
int toysBayesian; // number of toys used to average out nuisance paramereters in Bayesian method 
int seed;  // default 1234;   seed of random number generator engine
int systematics; // default 1, will determin if use systematics according to data card; if 0, then will not use systematics 
float CL;  // default 0.95;  Confidence Level...
lands::PRIOR prior; //default flat;   prior on signal.  other options:  corr,  1/sqrt(r)
int calcsignificance; // default 0;  if 1, then calc significance instead of limit
bool lowerLimit; // in FeldmanCousins,  0 by default for upper limit,  if 1 then calc lower limit  
int toysFactor; // in FeldmanCousins, Increase the toys per point by this factor w.r.t. the minimum from adaptive sampling,   default =1
bool adaptiveSampling; // currently only implemented in FeldmanCousins,  turn on (=1) by default.  =0 off
int testStat; //default LEP.  other options: TEV, Atlas, AtlasAllowMuHatNeg ...
int rule; //default CLs.  other options: CLsb
float clsAcc; // 0.001; Absolute accuracy on CLs to reach to terminate the scan
double rAbsAcc;// 0.01; Absolute accuracy on r to reach to terminate the scan
double rRelAcc;// 0.01; Relative accuracy on r to reach to terminate the scan

// for FeldmanCousins
bool bQuickEstimateInitialLimit = true; 
double initialRmin = 1., initialRmax = 21;// only when bQuickEstimateInitialLimit==false

int oneside = 1; //for PLR limit

int main(int argc, const char*argv[]){
	processParameters(argc, argv);

	CountingModel *cms; // this instance will be the one combining all sub cards
	cms = new CountingModel();
	cms->SetDebug(debug);

	/*
	 * combining at most 100 datacards
	 */
	CountingModel tmp1[100];
	CountingModel *tmp[100];// = new CountingModel(); 
	if(datacards.size()>100) {cout<<"too many datacards "<<datacards.size()<<endl; exit(0);}
	for(int i=0; i<datacards.size(); i++){
		tmp[i] = new CountingModel();
		tmp[i] -> SetDebug(debug);
		ConfigureModel(tmp[i], datacards[i].Data(), debug);
		tmp[i]->SetUseSystematicErrors(true);
	}
	if(debug)cout<<"totally "<<datacards.size()<<" data cards processed"<<endl;
	if(datacards.size()==1) cms=tmp[0];
	else if(datacards.size()>=2){
		tmp1[1] = CombineModels(tmp[0], tmp[1]);
		tmp1[1].SetUseSystematicErrors(true);
		if(debug)cout<<"2 data cards have been combined"<<endl;
		for(int i=2; i<datacards.size(); i++){
			tmp1[i] = CombineModels(&tmp1[i-1], tmp[i]);
			tmp1[i].SetUseSystematicErrors(true);
		}	
		cms = &tmp1[datacards.size()-1];
	}else{exit(0);}
	cout<<"totally "<<datacards.size()<<" data cards combined"<<endl;
	cms->SetUseSystematicErrors(systematics);
	// done combination

	// common operations
	if(debug)cms->Print();
	CRandom *rdm = new CRandom(seed);  //initilize a random generator
	cms->SetRdm(rdm);
	cms->SetUseSystematicErrors(systematics);
	cms->RemoveChannelsWithExpectedSignal0orBkg0();
	//cms->RemoveChannelsWithExpectedSignal0orBkg0(-1);
	//cms->SetAllowNegativeSignalStrength(false);
	if(dataset == "asimov_b")cms->UseAsimovData(0);
	else if(dataset == "asimov_sb")cms->UseAsimovData(1);


	// common results
	double rmean;
	vector<double> difflimits;

	TPaveText *pt = SetTPaveText(0.2, 0.7, 0.3, 0.9);
	if(calcsignificance==0){// calc limits
		if(method == "Bayesian"){
			BayesianBase bys(cms, 1-CL, 1.e-3);
			bys.SetDebug(debug);
			bys.SetNumToys(toysBayesian);
			bys.SetCrossSectionPrior(prior);
			double rtmp;
			rtmp = bys.Limit();
			cout<<"------------------------------------------------------------"<<endl;
			cout<<"Observed Upper Limit on the ratio R at 95\% CL = "<<rtmp<<endl;
			cout<<"------------------------------------------------------------"<<endl;

			// draw 

			if(debug)	{
				bys.PosteriorPdf();
				//start_time=cur_time; cur_time=clock(); cout << "\t TIME_in_BayesianLimit ppdf: "<< (cur_time - start_time)/1000. << " millisec\n";
				DrawEvolution2D pdfr(bys.GetVR(), bys.GetVP(), ";r;Likelihood", (jobname+"_postpdf").Data(), pt);	
				pdfr.setDrawPosteriorPdf(rtmp);
				pdfr.draw();
				pdfr.getGraph()->SetMarkerSize(0.01);
				pdfr.save();
			}

			if(doExpectation){
				cms->SetSignalScaleFactor(1.);
				LimitBands lb(&bys, cms);	
				lb.SetDebug(debug);
				int noutcomes = toys;
				lb.BysLimitBands(1-CL, noutcomes, toysBayesian);
				rmean=lb.GetBysLimitMean();
				difflimits = lb.GetDifferentialLimitsBys();

			}

		}else if(method == "Hybrid"){

			// initialize the calculator
			CLsBase frequentist;
			frequentist.SetDebug(debug);
			frequentist.SetRdm(rdm);
			frequentist.SetTestStatistics(testStat);
			cms_global= cms;
			//vdata_global=cms->Get_v_data();

			double tmp;
			vdata_global=cms->Get_v_data();


			frequentist.BuildM2lnQ(cms, toysHybrid);
			double errs, errb, errsb;
			double cls = frequentist.CLs(errs);
			double clsb = frequentist.CLsb(errsb);
			double clb = frequentist.CLb(errb);
			cout<<"------------------------------------------------------------"<<endl;
			cout<<"Observed CLs = "<<cls<<" +/- "<<errs<<endl;
			cout<<"Observed CLsb = "<<clsb<<" +/- "<<errsb<<endl;
			cout<<"Observed CLb = "<<clb<<" +/- "<<errb<<endl;
			cout<<"------------------------------------------------------------"<<endl;

			if(debug){
				DrawPdfM2logQ pdfM2logQ(frequentist.Get_m2logQ_sb(),frequentist.Get_m2logQ_b(), frequentist.Get_m2lnQ_data(), 
						"-2lnQ on data", "; -2lnQ; entries", "lnq", pt);
				pdfM2logQ.draw();
				double m2lnqdata = frequentist.Get_m2lnQ_data();
				//cout<<m2lnqdata<<endl;
				vector<double> vsb = frequentist.Get_m2logQ_sb();
				for(int i=0; i<vsb.size(); i++){
					//	if(vsb[i]>=m2lnqdata)cout<<i<<" "<<vsb[i]<<endl;
				}
			}

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
			clsr95.SetAlpha(1-CL);

			if(bQuickEstimateInitialLimit) rtmp = clsr95.LimitOnSignalScaleFactor(cms, &frequentist, toysHybrid);
			else rtmp = clsr95.LimitOnSignalScaleFactor(cms, initialRmin, initialRmax, &frequentist, toysHybrid, 3);

			if(debug){
				DrawEvolution2D d2d(clsr95.GetvTestedScaleFactors(), clsr95.GetvTestedCLs(), "; r ; CLs", (jobname+"_r_vs_cl").Data(), pt);
				d2d.draw();
				d2d.save();

			}

			cout<<"------------------------------------------------------------"<<endl;
			cout<<"Observed Upper Limit on the ratio R at 95\% CL = "<<rtmp<<" +/- "<<clsr95.LimitErr()<<endl;
			cout<<"------------------------------------------------------------"<<endl;

			if(doExpectation){
				cms->SetSignalScaleFactor(1.);
				LimitBands lb(&clsr95, &frequentist, cms);	
				lb.SetDebug(debug);
				int noutcomes = toys;
				lb.CLsLimitBands(1-CL, noutcomes, toysHybrid);
				rmean=lb.GetCLsLimitMean();
				difflimits = lb.GetDifferentialLimitsCLs();
			}

		}else if(method == "FeldmanCousins"){

			CLsBase frequentist;
			frequentist.SetDebug(debug);
			frequentist.SetRdm(rdm);
			cms_global= cms;
			//vdata_global=cms->Get_v_data();

			double tmp;
			vdata_global=cms->Get_v_data();

			// FeldmanCousins must use specified test statisics:  profile likelihood ratio, which allow mu hat > probed mu
			frequentist.SetTestStatistics(31);

			cms->SetAllowNegativeSignalStrength(false);


			CLsLimit clsr95;
			clsr95.SetDebug(debug);
			double rtmp;
			clsr95.SetAlpha(1-CL);


			clsr95.SetAdaptiveSampling(adaptiveSampling);
			bool lowerLimit_ = false;

			double r95_fc;

			double fcMid = 0, fcErr = 0; 
			double rmin = initialRmin,  rmax = initialRmax;

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
				r95_fc = clsr95.FeldmanCousins(cms, rmin, rmax, &frequentist, toysHybrid, nsteps);
				vv = clsr95.GetFCconstruction();
				for(int i=0; i<vv.size(); i++){
					vv_all.push_back(vv[i]);
				}
				int found = -1; 
				if (lowerLimit) {
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
				if (fcErr < 4*std::max(rAbsAcc, rRelAcc* fcMid)) { // make last scan more precise
					clsr95.SetAdditionalNToysFactor(4*toysFactor);
				}
			} while (fcErr > std::max(rAbsAcc, rRelAcc* fcMid));

			r95_fc = fcMid;
			if (0) {
				std::cout << "\n -- FeldmanCousins++ -- \n";
				std::cout << "Limit: r " << (lowerLimit_ ? "> " : "< ") << fcMid << " +/- "<<fcErr << "\n";
			}


			if(debug) cout <<"95\% CL upper limit by FC: "<<r95_fc<<",   use bys: "<<rtmp<<endl;

			cout<<"------------------------------------------------------------"<<endl;
			cout<<"Observed Upper Limit on the ratio R at 95\% CL = "<<r95_fc<<endl;
			cout<<"------------------------------------------------------------"<<endl;


			if(debug){
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
				delete r_tested;
				delete q_data;
			}

			return 1;
		}else if(method=="ProfiledLikelihood"){

			cms_global= cms;
			vdata_global=cms->Get_v_data();

			double r95;
			double tmp;
			double tmperr;
			double y0_1 =  MinuitFit(3, tmp, tmperr, 0);
			if(debug)	cout<<y0_1<<" fitter u="<<tmp<<" +/- "<<tmperr<<endl;
			double tmpr = 0;
			double y0_2 =  MinuitFit(2, tmpr, tmperr) ;
			if(debug)	cout<<y0_2<<" fitter u="<<tmp<<" +/- "<<tmperr<<endl;

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
			if (oneside==2) CI= 1.921;  // = 1.96**2/2 , two sided 95% CL --->  one sided 97.5%
			else CI = 1.64*1.64/2.; // two sided 90% CL
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

			if(debug)cout<<"r95 = "<<r95<<",  "<<nsearched<<" steps"<<endl;


			cout<<"------------------------------------------------------------"<<endl;
			cout<<"Observed Upper Limit on the ratio R at 95\% CL = "<<r95<<endl;
			cout<<"------------------------------------------------------------"<<endl;


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
				DrawEvolution2D profiledL(vrxsec, vmlnq, ";r = #sigma/#sigma_{SM}; -log#lambda(r)", (jobname+"_profiledL").Data(), pt);
				profiledL.setLogY(0);
				profiledL.draw();
			}

			return 1;

		}

		if(doExpectation){
			// plot the distribution of limits of noutcomes
			// and the cummulative pdf
			TString ts=jobname; ts+="_limits";
			FillTree(ts, difflimits);
			vector<double> all_calculated_R95s;
			vector<double> cummulativeProbabilities; // all_calculated_R95s has been sorted, each of them has a cummulative probability
			SortAndCumulative(difflimits, all_calculated_R95s, cummulativeProbabilities);

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
			printf("BANDS %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n", rm2s2, rm1s2, rmean, rp1s2, rp2s2, rmedian2);
			cout<<"------------------------------------------------------------"<<endl;

			TCanvas *c=new TCanvas("cme","cme");
			c->SetLogy(1);
			TH1F *h=new TH1F("h",";r=#frac{#sigma_{95%CL}}{#sigma_{SM}}; entries", 200, 0, 15);	
			for(int i=0; i<difflimits.size(); i++) h->Fill(difflimits[i]);
			h->Draw();
			Save(c, (jobname+"_differential_limits").Data());

			pt = SetTPaveText(0.67, 0.4, 0.8, 0.7);
			//pt->AddText("statistical bands");
			string ssave= (jobname+"_cummulativeR").Data();
			string stitle="; r = #sigma_{95%}/#sigma_{SM}; cumulative probability;";
			PlotXvsCummulativeProb plotRvsP(all_calculated_R95s, cummulativeProbabilities,
					rm1s2, rp1s2, rm2s2, rp2s2, ssave, stitle, pt);
			plotRvsP.draw();
			plotRvsP.getGraph()->SetMarkerSize(1);
			plotRvsP.save();

		}

	}else { // calc significances 

		// initialize the calculator
		CLsBase frequentist;
		frequentist.SetDebug(debug);
		frequentist.SetRdm(rdm);
		frequentist.SetTestStatistics(testStat);

		double tmp, tmperr;
		double rmean, rm1s, rm2s, rp1s, rp2s;	
		vector<double> difflimits; 
		if(method == "ProfiledLikelihood"){

			cms_global= cms;
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
			if(doExpectation){
				SignificanceBands lb(&frequentist, cms);	
				lb.SetDebug(debug);
				int noutcomes = toys;
				lb.Bands(noutcomes);
				rmean=lb.GetSignificanceMean();
				rm1s=lb.GetSignificance(-1);rm2s=lb.GetSignificance(-2);rp1s=lb.GetSignificance(1);rp2s=lb.GetSignificance(2);
				difflimits=lb.GetDifferentialSignificances();
			}
		}else if(method == "Hybrid"){

			cms->SetSignalScaleFactor(1.);
			int ntoysToDoSignificance = toysHybrid; //10000000;
			frequentist.SetModel(cms);
			cout<<"\n Running "<<ntoysToDoSignificance<<" toys to evaluate significance for data "<<endl;
			double signi = frequentist.SignificanceForData(ntoysToDoSignificance);
			if(debug){
				cout<<"Q_b_data = "<<frequentist.Get_m2lnQ_data()<<endl;	
				vector<double> vclb = frequentist.Get_m2logQ_b();
				TString  s = jobname; 
				s+="_hybridSig_ts"; s+=testStat;
				s+="_seed"; s+=seed;
				FillTree(s, vclb);

			}
			cout<<"------------------------------------------------------------"<<endl;
			cout<<" Observed Significance for the data = "<<signi<<endl;
			cout<<"------------------------------------------------------------"<<endl;

			if(doExpectation){
				// FIXME
			}
		}
		if(doExpectation){

			TString ts=jobname; ts+="_significances";
			FillTree(ts, difflimits);

			TCanvas *c=new TCanvas("csig","cSig");
			c->SetLogy(1);
			TH1F *h=new TH1F("h",";ProfiledLikelihood significance; entries", 200, 0, 8);	
			for(int i=0; i<difflimits.size(); i++) h->Fill(difflimits[i]);
			h->Draw();
			Save(c, (jobname+"differential_sig").Data());

			vector<double> all_calculated_R95s;
			vector<double> cummulativeProbabilities; // all_calculated_R95s has been sorted, each of them has a cummulative probability
			SortAndCumulative(difflimits, all_calculated_R95s, cummulativeProbabilities);

			pt = SetTPaveText(0.67, 0.4, 0.8, 0.7);
			pt->AddText("PL significance statistical bands");
			string ssave=(jobname+"_cummulativeS").Data();
			string stitle="; ProfiledLikelihood significance; cumulative probability;";
			PlotXvsCummulativeProb plotRvsP(all_calculated_R95s, cummulativeProbabilities,
					rm1s, rp1s, rm2s, rp2s,ssave, stitle, pt);
			plotRvsP.draw();

		}
	}
	return 1;
}

void processParameters(int argc, const char* argv[]){
	vector<TString> allargs;
	for(int i=1; i<argc; i++){
		allargs.push_back(TString(argv[i]));
	}

	for(int i=0; i<allargs.size(); i++){
		if(allargs[i]=="--help" or allargs[i]=="-h" ) PrintHelpMessage();
		//cout<<allargs[i]<<endl;
		if(allargs[i].BeginsWith("-") and !allargs[i].IsFloat()){
			TString key = allargs[i];
			vector<TString> values;
			for(int j=i+1; j<allargs.size(); j++){
				if(allargs[j].BeginsWith("-") and !allargs[j].IsFloat()) break;
				values.push_back(allargs[j]);
			}
			options.insert(TPair(key, values));
		}
	}

	TMap::iterator p;

	cout<<"Your arguments: "<<endl;
	for(p=options.begin(); p!=options.end();++p){
		cout<<p->first<<" ";
		for(int j=0; j<p->second.size(); j++){
			cout<<p->second[j]<<" ";
		}
		cout<<endl;
	}

	vector<TString> tmpcards = options["-d"];
	if(tmpcards.size()==0) tmpcards = options["--datacards"];
	cout<<endl<<"You are trying to combine following "<<tmpcards.size()<<" cards:"<<endl;
	if(tmpcards.size()==0) { cout<<"*** No data card specified, please use option \"-d or --datacards\""<<endl; exit(0); }
	else{
		for(int i=0; i<tmpcards.size(); i++){
			FileStat_t buf;
			if(gSystem->GetPathInfo(tmpcards[i], buf) ==1 ) { cout<<tmpcards[i] << " not found, skipped"<<endl;}
			else {
				cout<<tmpcards[i]<<endl;
				datacards.push_back(tmpcards[i]);
			}
		}
		cout<<" .... valid data cards = "<<datacards.size()<<endl;
		if(datacards.size()<=0){ cout<< " please provide valid data cards "<<endl; exit(0); }
	}

	vector<TString> tmpv;
	// limit or significance
	tmpv = options["--significance"]; 
	if( tmpv.size()!=1 ) { calcsignificance = 0; }
	else calcsignificance = tmpv[0].Atoi();

	tmpv = options["-M"]; if(tmpv.size()!=1) tmpv = options["--method"];
	if( tmpv.size()!=1 ) { cout<<"ERROR No method specified, please use option \"-M or --method\" "<<endl; exit(0); }
	else {
		method = tmpv[0];
		if( calcsignificance && 
				(method!="ProfiledLikelihood" and method!="Hybrid")
		  ){cout<<"ERROR You are trying to use "<<method<<", which is not supported currently to calculate Significance"<<endl; exit(0);}
		if( !calcsignificance && 
				(method!="ProfiledLikelihood" and method!="Hybrid" and method!="Bayesian" and method!="FeldmanCousins" )
		  ){cout<<"ERROR You are trying to use "<<method<<", which is not supported currently to calculate limit "<<endl; exit(0);}
	}

	tmpv = options["-v"]; if(tmpv.size()!=1) tmpv = options["--verbose"]; if(tmpv.size()!=1) tmpv = options["--debug"]; 
	if( tmpv.size()!=1 ) { debug = 0; }
	else debug = tmpv[0].Atoi();

	tmpv = options["-n"]; if(tmpv.size()!=1) tmpv = options["--name"];
	if( tmpv.size()!=1 ) { jobname = datacards[0]+"_"+method; }
	else jobname = tmpv[0];

	tmpv = options["-D"]; if(tmpv.size()!=1) tmpv = options["--dataset"];
	if( tmpv.size()!=1 ) { dataset = "data_obs"; }
	else { 
		dataset = tmpv[0];
		if(dataset!="data_obs" and dataset!="asimov_sb" and dataset!="asimov_b"){cout<<"ERROR: dataset option must be one of data_obs, asimov_sb and asimov_b"<<endl; exit(0);}
	}

	tmpv = options["--doExpectation"]; 
	if( tmpv.size()!=1 ) { doExpectation = 0; }
	else doExpectation = tmpv[0].Atoi();

	tmpv = options["-t"]; if(tmpv.size()!=1) tmpv = options["--toys"];
	if( tmpv.size()!=1 ) { toys = 1000; }
	else toys = tmpv[0].Atoi();

	tmpv = options["-s"]; if(tmpv.size()!=1) tmpv = options["--seed"];
	if( tmpv.size()!=1 ) { seed = 1234; }
	else seed = tmpv[0].Atoi();

	tmpv = options["-S"]; if(tmpv.size()!=1) tmpv = options["--systematics"];
	if( tmpv.size()!=1 ) { systematics = 1; }
	else systematics = tmpv[0].Atoi();

	tmpv = options["-C"]; if(tmpv.size()!=1) tmpv = options["--CL"];
	if( tmpv.size()!=1 ) { CL = 0.95; }
	else {
		CL = tmpv[0].Atof();
		if(CL<=0 or CL>=1) {cout<<"ERROR  CL must be in the range (0, 1)"<<endl; exit(0);}
	}

	// Bayesian specific options
	tmpv = options["--prior"];
	if( tmpv.size()!=1 ) { prior = flat; }
	else {
		if(tmpv[0]=="flat"){prior=flat;}
		else if(tmpv[0]=="corr"){prior=corr;}
		else if(tmpv[0]=="1/sqrt(r)"){prior=prior_1overSqrtS;}
		else { cout<<"ERROR  supported baysian signal priors: flat , corr, 1/sqrt(r) "<<endl; exit(0);}
	}

	tmpv = options["-tB"]; if(tmpv.size()!=1) tmpv = options["--toysBayesian"];
	if( tmpv.size()!=1 ) { toysBayesian = 10000; }
	else {
		toysBayesian = tmpv[0].Atoi();
		if(toysBayesian<0) toysBayesian = 1;
	}


	// FeldmanCousins specific options
	tmpv = options["--lowerLimit"]; 
	if( tmpv.size()!=1 ) { lowerLimit = false; }
	else {
		if(tmpv[0].Atoi()==0) lowerLimit=false;
		else lowerLimit = true;
	}

	tmpv = options["--toysFactor"]; 
	if( tmpv.size()!=1 ) { toysFactor = 1; }
	else {
		toysFactor = tmpv[0].Atoi();
		if(toysFactor<1) toysFactor = 1;
	}

	tmpv = options["--adaptiveSampling"]; 
	if( tmpv.size()!=1 ) { adaptiveSampling = true; }
	else {
		if(tmpv[0].Atoi()==0) adaptiveSampling=false;
		else adaptiveSampling = true;
	}

	tmpv = options["--bQuickEstimateInitialLimit"]; 
	if( tmpv.size()!=1 ) { bQuickEstimateInitialLimit = true; }
	else {
		if(tmpv[0].Atoi()==0) bQuickEstimateInitialLimit=false;
		else bQuickEstimateInitialLimit = true;
	}

	tmpv = options["--initialRmin"]; 
	if( tmpv.size()!=1 ) { initialRmin = 1; }
	else {
		initialRmin = tmpv[0].Atof();
		//if(initialRmin<0) initialRmin = 1;
	}

	tmpv = options["--initialRmax"]; 
	if( tmpv.size()!=1 ) { initialRmax = 20; }
	else {
		initialRmax = tmpv[0].Atof();
		//if(initialRmax<=0) initialRmax = 20;
	}


	// Hybrid specific options
	tmpv = options["-tH"]; if(tmpv.size()!=1) tmpv = options["--toysHybrid"];
	if( tmpv.size()!=1 ) { toysHybrid = 10000; }
	else {
		toysHybrid = tmpv[0].Atoi();
		if(toysHybrid<0) toysHybrid = 100;
	}

	tmpv = options["--testStat"]; 
	if( tmpv.size()!=1 ) { testStat = 1; }
	else {
		if(tmpv[0]=="LEP") testStat=1;
		else if(tmpv[0]=="TEV") testStat=2;
		else if(tmpv[0]=="Atlas") testStat=3;
		else if(tmpv[0]=="AtlasAllowMuHatNeg") testStat=32;
		else {cout<<"ERROR Unimplemented testStat: "<<tmpv[0]<<". Supported: LEP, TEV, Atlas, AtlasAllowMuHatNeg "<<endl; exit(0); }
	}

	tmpv = options["--rule"]; 
	if( tmpv.size()!=1 ) { rule = 1; }
	else { 
		if(tmpv[0]=="CLs") rule=1;
		else if(tmpv[0]=="CLsb" or tmpv[0]=="CLsplusb") rule=2;
		else {cout<<"ERROR Unimplemented rule: "<<tmpv[0]<<". Supported: CLs,  CLsb "<<endl; exit(0); }
	}

	tmpv = options["--clsAcc"];
	if( tmpv.size()!=1 ) { clsAcc = 0.001; }
	else clsAcc = tmpv[0].Atof();

	tmpv = options["--rAbsAcc"];
	if( tmpv.size()!=1 ) { rAbsAcc = 0.01; }
	else rAbsAcc = tmpv[0].Atof();

	tmpv = options["--rRelAcc"];
	if( tmpv.size()!=1 ) { rRelAcc = 0.01; }
	else rRelAcc = tmpv[0].Atof();

	// ProfiledLikelihood specific options
	tmpv = options["--OneOrTwoSided"];
	if( tmpv.size()!=1 ) { oneside = 1; }
	else {
		oneside = tmpv[0].Atoi();
		if(oneside!=1 and oneside!=2)  {cout<<"ERROR --OneOrTwoSided can only have 1 or 2 as arg"<<endl; exit(0);}
	}



	printf("\n\n[ Summary of configuration in this job: ]\n");
	cout<<"  Calculating "<<(calcsignificance?"significance":"limit")<<" with "<<method<<" method "<<endl;
	cout<<"  datacards: "; for(int i=0; i<datacards.size(); i++) cout<<datacards[i]<<" "; cout<<endl;
	cout<<"  "<<(systematics?"use systematics":"not use systematics")<<endl;
	cout<<"  dataset is "<<dataset<<endl;
	if(!calcsignificance) cout<<"  target confidence level = "<<CL<<endl;
	cout<<"  "<<(doExpectation?"also calc expectation bands":"do not calc expectation bands")<<endl;
	if(doExpectation)cout<<"  number of outcomes to build expecation bands: "<<toys<<endl;
	if(method=="Bayesian") { 
		cout<<"  prior = ";
		if(prior==flat) cout<<"flat"<<endl;
		else if(prior==corr) cout<<"corr"<<endl;
		else if(prior==prior_1overSqrtS) cout<<"1/sqrt(r)"<<endl;
		else cout<<"Unknow"<<endl;

		cout<<"  toysBayesian = "<<toysBayesian<<endl;
	}else if(method=="ProfiledLikelihood"){
		if(!calcsignificance)cout<<(oneside==1?"  one sided":"  two sided")<<endl;
	}else if(method=="FeldmanCousins"){
		cout<<"  calc "<<(lowerLimit?"lower limit":"upper limit")<<endl;
		cout<<"  "<<(adaptiveSampling?"use adaptiveSampling":"do not use adaptiveSampling")<<endl;
		cout<<"  toysFactor = "<<toysFactor<<endl;
		if(adaptiveSampling==false) cout<<"   number of toys per point = "<<toysHybrid<<endl;
		cout<<"  do "<<(bQuickEstimateInitialLimit?"":"not")<<" estimate initial limit with bayesian method"<<endl;
		if(!bQuickEstimateInitialLimit){
			cout<<"   initialRmin = "<<initialRmin<<", initialRmax = "<<initialRmax<<endl;
			if(initialRmax<initialRmin) {cout<<"  initialRmin can't be > initialRmax"<<endl; exit(0);}
		}
	}else if(method=="Hybrid"){
		if(calcsignificance==false){
			cout<<"  testStat = "; if(testStat==1) cout<<"LEP"; if(testStat==2)cout<<"TEV"; if(testStat==3)cout<<"Atals"; if(testStat==32)cout<<"AtlasAllowMuHatNeg";
			cout<<endl;
			cout<<"  rule     = "; if(rule==1) cout<<"CLs"; if(rule==2)cout<<"CLsb";
			cout<<endl;
		}
		cout<<"  number of toys to build Q distribution = "<<toysHybrid<<endl;
	}

	cout<<"  random number generator seed: "<<seed<<endl;
	cout<<"  debug level = "<<debug<<endl;
	cout<<"  job name: "<<jobname<<endl;
	cout<<endl<<endl;

	fflush(stdout);	
}

void PrintHelpMessage(){
	printf("Usage: ./lands.exe [options] \n");                                                                                                                 
	printf("Allowed options: \n");                                                                                                                 
	printf("-h [ --help ]                         Produce help message \n"); 
	printf("-v [ --verbose ] arg (=0)             Verbosity level \n"); 
	printf("-n [ --name ] arg                     Name of the job,  default is \"datacard\"+\"method\" \n"); 
	printf("-d [ --datacards ] args               Datacard files,  can contain \"*, ?\" \n"); 
	printf("-D [ --dataset ] arg (=data_obs)      Dataset for observed limit,  data_obs,  asimov_b, asimov_sb \n"); 
	printf("-M [ --method ] arg                   Method to extract upper limit. Supported methods are: Bayesian, FeldmanCousins, Hybrid, ProfiledLikelihood \n"); 
	printf("--doExpectation arg (=0)              i.e calc expected bands and mean/median values     \n"); 
	printf("-t [ --toys ] arg (=1000)             Number of Toy MC extractions for expectation \n"); 
	printf("-s [ --seed ] arg (=1234)             Toy MC random seed \n"); 
	printf("-S [ --systematics ] arg (=1)         if 0, then will not use systematics  \n"); 
	printf("-C [ --cl ] arg (=0.95)               Confidence Level \n"); 
	printf("--significance arg (=0)               Compute significance instead of upper limit,  supported methods: ProfiledLikelihood and Hybrid (CLb) \n"); 
	printf(" \n"); 
	printf("Bayesian specific options: \n"); 
	printf("--prior arg (=flat)            	      Prior to use: \'flat\' (default), \'1/sqrt(r)\', \'corr\' \n"); 
	printf("-tB [ --toysBayesian ] arg (=10000)   number of toys used to average out nuisance paramereters in Bayesian method     \n"); 
	printf(" \n"); 
	printf("FeldmanCousins specific options: \n"); 
	printf("--lowerLimit arg (=0)                 Compute the lower limit instead of the upper limit \n"); 
	printf("--toysFactor arg (=1)                 Increase the toys per point by this factor w.r.t. the minimum from adaptive sampling \n"); 
	printf("--adaptiveSampling arg (=1)           currently only implemented in FeldmanCousins,  turn on (=1) by default.  =0 off \n"); 
	printf("--bQuickEstimateInitialLimit arg (=1) quickly estimate initial limit from bayesian technique, turn off by 0  \n"); 
	printf("--initialRmin arg (=1)                only effective when bQuickEstimateInitialLimit=0 \n"); 
	printf("--initialRmax arg (=20)               only effective when bQuickEstimateInitialLimit=0 \n"); 
	printf(" \n"); 
	printf("Hybrid specific options: \n"); 
	printf("-tH [ --toysHybrid ] arg (=10000)     Number of Toy MC extractions to compute CLs+b, CLb and CLs \n"); 
	printf("--clsAcc arg (=0.001)                 Absolute accuracy on CLs to reach to terminate the scan \n"); 
	printf("--rAbsAcc arg (=0.01)                 Absolute accuracy on r to reach to terminate the scan \n"); 
	printf("--rRelAcc arg (=0.01)                 Relative accuracy on r to reach to terminate the scan \n"); 
	printf("--rule arg (=CLs)                     Rule to use: CLs, CLsb \n"); 
	printf("--testStat arg (=LEP)                 Test statistics: LEP, TEV, Atlas, AtlasAllowMuHatNeg. \n"); 
	printf(" \n"); 
	printf("ProfiledLikelihood specific options: \n"); 
	printf("--OneOrTwoSided arg (=1)              1 sided limit -lnL = 1.345;  2 sided limit -lnL = 1.921 \n"); 
	exit(0);
}
