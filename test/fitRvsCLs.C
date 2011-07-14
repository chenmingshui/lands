/*
 * =====================================================================================
 *
 *       Filename:  fitRvsCLs.C
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/01/2011 10:23:18 AM CEST
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Mingshui Chen (), Mingshui.Chen@cern.ch
 *        Company:  Univsity of Florida
 *
 * =====================================================================================
 */
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TFile.h"
#include "TObject.h"
#include "TTree.h"
#include "TH1D.h"
#include "TMath.h"
#include "TPaveText.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TKey.h"
#include "TCanvas.h"

#include <iostream>
#include <string>
#include <map>
#include <utility>
#include <vector>
#include <iomanip>
#include <math.h>
#include <ctime> // upto second
#include <time.h> // upto micro second
#include <stdexcept>
#include <algorithm>


void readAllToysFromFile(TGraphErrors* tge, TFile*f, double d_quantile= -1); // when quantile <=0: make distribution for observed data.   
bool GetCLs(double qdata, TTree* tsb, TTree*tb,  double &cls, double &err, double d_quantile=-1);
bool GetM2lnQ(TTree* tsb, TTree*tb, vector<double> &vqsb, vector<double>&vqb);
bool GetPValue(vector<double> vqsb, double qdata, double &ret, double &err);
int _debug = 0;
double rMin, rMax;
typedef std::pair<double,double> CLs_t;
CLs_t clsMin, clsMax;
double clsTarget = 0.05;
double limit, limitErr;
double rAbsAccuracy_ = 0.01;
double rRelAccuracy_ = 0.02;
bool bPlotInFittedRange = false;
double extractLimitAtQuantile(TString inFileName, TString plotName, double d_quantile = -1 );
void fitRvsCLs(){
}
void run(TString inFileName, TString plotName){
	double m2s =	extractLimitAtQuantile(inFileName, plotName+"_-2sigma", 0.0275 );
	double m1s =	extractLimitAtQuantile(inFileName, plotName+"_-1sigma", 0.16 );
	double med =	extractLimitAtQuantile(inFileName, plotName+"_median", 0.5 );
	double p1s =	extractLimitAtQuantile(inFileName, plotName+"_1sigma", 0.84 );
	double p2s =	extractLimitAtQuantile(inFileName, plotName+"_2sigma", 0.975 );
	double dat =	extractLimitAtQuantile(inFileName, plotName+"_observed", -1 );

	cout<<"EXPECTED LIMIT BANDS from(-2s,-1s,median,1s,2s): "<<m2s<<" "<<m1s<<" "<<med<<" "<<p1s<<" "<<p2s<<endl;
	cout<<"Observed data limit: "<<dat<<endl;
}
double extractLimitAtQuantile(TString inFileName, TString plotName, double d_quantile ){
	TFile *f = TFile::Open(inFileName);
	TF1 *expoFit = new TF1("expoFit","[0]*exp([1]*(x-[2]))", rMin, rMax);
	TGraphErrors *limitPlot_ =  new TGraphErrors();

	bool done = false;
	if (_debug > 0) std::cout << "Search for upper limit using pre-computed grid of p-values" << std::endl;

	readAllToysFromFile(limitPlot_, f, d_quantile ); 
	limitPlot_->Sort();
	double minDist=1e3;
	int n= limitPlot_->GetN();
	cout<<" Number of points in limitPlot_ : "<<n<<endl;

	clsMin.first=0;
	clsMin.second=0;
	clsMax.first=0;
	clsMax.second=0;

	for (int i = 0; i < n; ++i) {
		double x = limitPlot_->GetX()[i], y = limitPlot_->GetY()[i], ey = limitPlot_->GetErrorY(i);
		if (_debug > 0) std::cout << "  r " << x << ", CLs = " << y << " +/- " << ey << std::endl;
		if (y-3*ey >= clsTarget) { rMin = x; clsMin = CLs_t(y,ey); }
		if (y+3*ey <= clsTarget) { rMax = x; clsMax = CLs_t(y,ey); }
		if (fabs(y-clsTarget) < minDist) { limit = x; minDist = fabs(y-clsTarget); }
	}
	if((clsMin.first==0 and clsMin.second==0) || (clsMax.first==0 and clsMax.second==0))
	{
		rMin = limitPlot_->GetX()[0]; clsMin=CLs_t(limitPlot_->GetY()[0], limitPlot_->GetErrorY(0)); 
		rMax = limitPlot_->GetX()[n-1]; clsMax=CLs_t(limitPlot_->GetY()[n-1], limitPlot_->GetErrorY(n-1));
	}

	if (_debug > 0) std::cout << " after scan x ~ " << limit << ", bounds [ " << rMin << ", " << rMax << "]" << std::endl;
	limitErr = std::max(limit-rMin, rMax-limit);
	expoFit->SetRange(rMin,rMax);
	//expoFit.SetRange(limitPlot_->GetXaxis()->GetXmin(),limitPlot_->GetXaxis()->GetXmax());
	//expoFit->SetRange(1.7,2.25);

	if (limitErr < std::max(rAbsAccuracy_, rRelAccuracy_ * limit)) {
		if (_debug > 1) std::cout << "  reached accuracy " << limitErr << " below " << std::max(rAbsAccuracy_, rRelAccuracy_ * limit) << std::endl;
		done = true; 
	}

	//if (!done) { // didn't reach accuracy with scan, now do fit
	if (1) { // didn't reach accuracy with scan, now do fit
		if (_debug) {
			std::cout << "\n -- HybridNew, before fit -- \n";
			std::cout << "Limit: r" << " < " << limit << " +/- " << limitErr << " [" << rMin << ", " << rMax << "]\n";

			std::cout<<"rMin="<<rMin<<" clsMin="<<clsMin.first<<",  rMax="<<rMax<<" clsMax="<<clsMax.first<<endl;
		}

		expoFit->FixParameter(0,clsTarget);
		expoFit->SetParameter(1,log(clsMax.first/clsMin.first)/(rMax-rMin));
		expoFit->SetParameter(2,limit);
		double rMinBound, rMaxBound; expoFit->GetRange(rMinBound, rMaxBound);
		limitErr = std::max(fabs(rMinBound-limit), fabs(rMaxBound-limit));
		int npoints = 0; 
		for (int j = 0; j < limitPlot_->GetN(); ++j) { 
			if (limitPlot_->GetX()[j] >= rMinBound && limitPlot_->GetX()[j] <= rMaxBound) npoints++; 
		}
		for (int i = 0, imax = 0; i <= imax; ++i, ++npoints) {
			limitPlot_->Sort();
			limitPlot_->Fit(expoFit,(_debug <= 1 ? "QNR EX0" : "NR EXO"));
			if (_debug) {
				std::cout << "Fit to " << npoints << " points: " << expoFit->GetParameter(2) << " +/- " << expoFit->GetParError(2) << std::endl;
			}

			// only when both  "cls+3e<0.05 and cls-3e>0.05" are satisfied, we require below ...
			//	if ((rMin < expoFit->GetParameter(2))  && (expoFit->GetParameter(2) < rMax) && (expoFit->GetParError(2) < 0.5*(rMaxBound-rMinBound))) { 
			// sanity check fit result
			limit = expoFit->GetParameter(2);
			limitErr = expoFit->GetParError(2);
			if (limitErr < std::max(rAbsAccuracy_, rRelAccuracy_ * limit)) break;
			//	}
		} 
	}

	if (limitPlot_) {
		TCanvas *c1 = new TCanvas("c1","c1");
		limitPlot_->Sort();
		limitPlot_->SetLineWidth(2);
		/*
		   double xmin = 0, xmax = 10; //FIXME
		   for (int j = 0; j < limitPlot_->GetN(); ++j) {
		   if (limitPlot_->GetY()[j] > 1.4*clsTarget || limitPlot_->GetY()[j] < 0.6*clsTarget) continue;
		   xmin = std::min(limitPlot_->GetX()[j], xmin);
		   xmax = std::max(limitPlot_->GetX()[j], xmax);
		   }
		   */
		double rMinBound, rMaxBound; expoFit->GetRange(rMinBound, rMaxBound);
		if(bPlotInFittedRange){
			limitPlot_->GetXaxis()->SetRangeUser(rMinBound, rMaxBound);
			limitPlot_->GetYaxis()->SetRangeUser(0.5*clsTarget, 1.5*clsTarget);
		}
		limitPlot_->Draw("AP");
		expoFit->Draw("SAME");
		TLine line(limitPlot_->GetX()[0], clsTarget, limitPlot_->GetX()[limitPlot_->GetN()-1], clsTarget);
		line.SetLineColor(kRed); line.SetLineWidth(2); line.Draw();
		line.DrawLine(limit, 0, limit, limitPlot_->GetY()[0]);
		line.SetLineWidth(1); line.SetLineStyle(2);
		line.DrawLine(limit-limitErr, 0, limit-limitErr, limitPlot_->GetY()[0]);
		line.DrawLine(limit+limitErr, 0, limit+limitErr, limitPlot_->GetY()[0]);
		c1->Print(plotName+".gif");
	}

	std::cout << "\n -- Hybrid New -- \n";
	std::cout << "Limit: r" << " < " << limit << " +/- " << limitErr << " @ " << (1-clsTarget) * 100 << "% CL\n";
	return limit;
}
void readAllToysFromFile(TGraphErrors*tge, TFile*f, double d_quantile) {
	tge->Set(0);	
	if (!f) throw std::logic_error("Cannot use readHypoTestResult: option toysFile not specified, or input file empty");
	TDirectory *toyDir = f->GetDirectory("");
	if (!toyDir) throw std::logic_error("Cannot use readHypoTestResult: empty toy dir in input file empty");
	TIter next(toyDir->GetListOfKeys()); TKey *k;
	std::map<double, TTree*> gridCLsb; //r, <clsb, clserr>
	std::map<double, TTree*> gridCLb; //r, <clsb, clserr>
	std::map<double, double> gridQdata; //r, q_data
	while ((k = (TKey *) next()) != 0) {
		double rVal;
		TString name(k->GetName());
		if(name.BeginsWith("SAMPLING_SB_TESTED_R")){
			rVal = name.ReplaceAll("SAMPLING_SB_TESTED_R","").Atof();
			TTree *t = dynamic_cast<TTree *>(toyDir->Get(k->GetName()));
			gridCLsb[rVal]=t;
			if (_debug > 2) std::cout << "  Do " << k->GetName() << " -> " << name << " --> " << rVal << " tsb->entries= "<<t->GetEntries()<< std::endl;
		}else if(name.BeginsWith("SAMPLING_B_TESTED_R")){
			rVal = name.ReplaceAll("SAMPLING_B_TESTED_R","").Atof();
			TTree *t = dynamic_cast<TTree *>(toyDir->Get(k->GetName()));
			gridCLb[rVal]=t;
			if (_debug > 2) std::cout << "  Do " << k->GetName() << " -> " << name << " --> " << rVal << " tb->entries= "<<t->GetEntries()<< std::endl;
		}else if(name.BeginsWith("DATA_R")){
			name.ReplaceAll("DATA_R","");
			TString tmp = name;
			rVal = tmp.Remove(tmp.Index("_"), tmp.Length()).Atof();
			name.Remove(0,name.Index("_Q")+2);
			gridQdata[rVal]=name.Atof();
			if (_debug > 2) std::cout << "  Do " << k->GetName() << " -> " << tmp << " --> " << rVal << " Q_data ="<<name<<" --> "<<name.Atof()<< std::endl;
		}

	}

	int i = 0, n = gridCLsb.size();
	cout<<" grid size = "<<n<<endl;
	int n_valid = 0;
	for (std::map<double, TTree *>::iterator itg = gridCLsb.begin(), edg = gridCLsb.end(); itg != edg; ++itg) {
		double cls, clserr;
		GetCLs(gridQdata[itg->first], itg->second, gridCLb[itg->first], cls, clserr, d_quantile);
		//if(itg->first != 1.74 && itg->first != 1.95) continue;
		if((cls!=1 and clserr==0) or cls<1e-10)continue;
		n_valid+=1;
		tge->Set(n_valid);
		tge->SetPoint(     n_valid-1, itg->first, cls   ); 
		tge->SetPointError(n_valid-1, 0,          clserr);
		cout<<" input grid:  r="<<itg->first<<" cls="<<cls<<"+/-"<<clserr<<endl;
		i++;
	}
	cout<<"tge->N = "<<tge->GetN()<<endl;
}
bool GetCLs(double qdata, TTree* tsb, TTree*tb,  double &cls, double &err, double d_quantile){
	vector<double> vqsb, vqb; vqsb.clear(); vqb.clear();
	GetM2lnQ(tsb, tb, vqsb, vqb);

	if(d_quantile>0 && d_quantile<1){
		std::sort(vqb.begin(), vqb.end());
		int pos = int(vqb.size()*(1-d_quantile));// in -2lnQ distribution, the more left, the more signal like
		qdata = vqb[(pos<vqb.size()?pos:(vqb.size()-1))];
	}

	double clsb, clb, errsb, errb;
	GetPValue(vqsb, qdata, clsb, errsb);
	GetPValue(vqb, qdata, clb, errb);

	if(_debug)cout<<"qdata = "<<qdata<<endl;
	if(clb==0){if(_debug)	cout<<"CLsBase::CLs  Warning clb_b==0 !!!!"<<endl; err = 1;  return 1;}
	err = sqrt( errb/clb*errb/clb + errsb/clsb*errsb/clsb) * clsb/clb;
	if(_debug) cout<<"CLsb = "<<clsb<<"+/-"<<errsb<<endl;
	if(_debug) cout<<"CLb = "<<clb<<"+/-"<<errb<<endl;
	if(_debug) cout<<"CLsBase::CLs  CLs=CLsb/CLb="<<clsb/clb<<"+/-"<<err<<endl;
	cls = clsb/clb;

	return true;
}
bool GetM2lnQ(TTree* tsb, TTree*tb, vector<double> &vqsb, vector<double>&vqb){
	double clsb, clb;
	TBranch *brCLsb ;
	tsb->SetBranchAddress("brT", &clsb, &brCLsb);
	Long64_t nentries = tsb->GetEntries();
	for(int i=0; i<nentries; i++){
		tsb->GetEntry(i);
		vqsb.push_back(clsb);
	}
	//for(int i=0; i<vqsb.size(); i++) cout<<" "<<vqsb[i]<<" "<<endl;
	TBranch *brCLb ;
	tb->SetBranchAddress("brT", &clb, &brCLb);
	nentries = tb->GetEntries();
	for(int i=0; i<nentries; i++){
		tb->GetEntry(i);
		vqb.push_back(clb);
	}
	//for(int i=0; i<vqb.size(); i++) cout<<" "<<vqb[i]<<" "<<endl;
	return true;		
}
bool GetPValue(vector<double> vqsb, double qdata, double &ret, double &err){
	double tmp = qdata;
	int _nexps = vqsb.size();
	for(int i=0; i<_nexps; i++){
		if(vqsb[i]>=tmp)
			ret ++ ;	
	}		
	ret/=_nexps;

	err= sqrt(ret*(1-ret)/_nexps);
	if(ret==0||ret==1) err= 1./_nexps;

	if(_debug>=10){
		cout<<"p ="<<ret<<" +/- "<<err<<" and total exps="<<_nexps<<endl;
		if(ret*_nexps <= 20) cout<<"p*nexps="<<ret*_nexps<<", statistic may not enough"<<endl;
	}
	if(ret == 0){
		if(_debug)	cout<<"p=0, it means number of pseudo experiments is not enough"<<endl;
		if(_debug)	cout<<"              Currently, we put p=1./"<<_nexps<<endl;
		ret = 1./(double)_nexps;
	}
	return true;
}
