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
//#include "inputs.C"

#include "TGraph.h"

#include <time.h> // upto micro second

//#include "TPython.h"
//
#include "TString.h"
#include "TObjString.h"
#include "TObject.h"
#include "TObjArray.h"
#include "TIterator.h"
#include <fstream>
#include "TMinuit.h"
#include "TSystem.h"

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
double L_fit(int model, double &r , double &er  ){
	bool debugMinuit = 0;
	bool UseMinos = 0;

	int npars = cms_global->Get_max_uncorrelation();
	TMinuit *gMinuit = new TMinuit(npars+2);  //initialize TMinuit with a maximum of 5 params
	gMinuit->SetFCN(lands::Chisquare);

	Double_t arglist[10];
	Int_t ierflg = 0;


	if(!debugMinuit){
		arglist[0]=-1;
		gMinuit -> mnexcm("SET PRINT", arglist, 1, ierflg);
		gMinuit -> mnexcm("SET NOW", arglist, 1, ierflg);
		
	}

	arglist[0] = 2;
	//gMinuit->mnexcm("SET STRATEGY", arglist ,1,ierflg);
	//gMinuit -> mnexcm("SET NOG", arglist, 1, ierflg);


	// Set starting values and step sizes for parameters
	// gMinuit->mnparm(par_index, "par_name", start_value, step_size, lower, higher, ierflg);
	vector<int> v_pdftype = cms_global->Get_v_pdftype();
	vector<double> v_TG_maxUnc = cms_global->Get_v_TruncatedGaussian_maxUnc();
	for(int i=1; i<=cms_global->Get_max_uncorrelation(); i++){
		TString sname; 
		sname.Form("p%d",i);
		if(v_pdftype[i] == typeLogNormal )
			gMinuit->mnparm(i, sname, 0, 0.1, -5, 5,ierflg);
		else if(v_pdftype[i] == typeTruncatedGaussian ){
			double maxunc = v_TG_maxUnc[i];	
			if(maxunc>0.2) maxunc = -1./maxunc;
			else maxunc = -5;
			gMinuit->mnparm(i, sname, 0, 0.1, maxunc, 5,ierflg);
		}else {
			cout<<"pdftype not yet defined "<<endl;
			exit(0);
		}
	}

	// through fixing the ratio to determine whether fit for S+B(r=1) or B-only (r=0)   Q_tevatron
	// let the ratio float, then it's Q_atlas
	if(model==1){ // S+B, fix r
		gMinuit->mnparm(0, "ratio", 1, 0.1, 0, 100, ierflg);
		gMinuit->FixParameter(0);
	}
	else if(model==0){ // B-only, fix r
		gMinuit->mnparm(0, "ratio", 0.0, 0.1, -1, 100, ierflg);
		gMinuit->FixParameter(0);
	}
	else if(model==2){ // S+B,  float r
		gMinuit->mnparm(0, "ratio", 1, 0.1, 0, 100, ierflg);
	}
	else if(model==3){ // B,  float r
		gMinuit->mnparm(0, "ratio", 0.0, 0.1, -1, 100, ierflg);
		gMinuit->FixParameter(0);
	}else {
		cout<<"Model not specified correctly:  0-3"<<endl;
		return 0;
	}

	arglist[0] = 1;
	gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
	// Now ready for minimization step
	arglist[0] = 500;
	arglist[1] = 1.;
	gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
	//	gMinuit->mnexcm("MINI", arglist ,2,ierflg);
	//	gMinuit->mnexcm("IMPROVE", arglist ,2,ierflg);


	if(UseMinos){
		arglist[0] = 500;
		gMinuit->mnexcm("MINOS", arglist ,1,ierflg);

	}

	// Print results
	Double_t amin,edm,errdef;
	Int_t nvpar,nparx,icstat;
	gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
	//cout<<"Minimized L = "<<gMinuit->fAmin<<endl;
	double l = gMinuit->fAmin;
	gMinuit->GetParameter(0, r, er);
	delete gMinuit;
	return l;
}	
int main(int argc, const char* argv[]){

	//gSystem->Load("libMinuit");
	CRandom *rdm = new CRandom(1234);
	CountingModel* cms=new CountingModel();
	cms->SetRdm(rdm);

	/*
	   cms->AddChannel(1.2, 0.9);
	   cms->AddObservedData(0, 2);
	   cms->AddUncertainty(0, 0, 0.2, 1, 1);
	   cms->AddUncertainty(0, 1, 0.5, 1, 2);
	   */

	const char* fileName =  "/data/Projects/Statistics/LandS/test/ATLAS_CMS/cmsinput/2010.06.22_hww_mh140A_3ch_0jets.txt";
	if(argc>=2) fileName = argv[1];
	ConfigureModel(cms, fileName); 


	cms->SetUseSystematicErrors(true);
	cms->Print();


	clock_t start_time=clock(), cur_time=clock();
	cms_global= cms;
	vdata_global = cms->Get_v_data();
	double r_sb, er_sb; 
	double L_sb = L_fit(1, r_sb, er_sb);
	double r_b, er_b; 
	double L_b = L_fit(0, r_b, er_b);
	double r_sb_float, er_sb_float; 
	double L_sb_float = L_fit(2, r_sb_float, er_sb_float);
	double del;
	double L_b_float = L_fit(3, del, del );
	L_b_float = L_fit(3, del, del );
	L_b_float = L_fit(3, del, del );
	L_b_float = L_fit(3, del, del );
	L_b_float = L_fit(3, del, del );
	start_time=cur_time; cur_time=clock(); cout << "\t\t\t  " << (cur_time - start_time)/1000. << " millisec\n"; 
	cout<<"L_sb = "<<L_sb<< ",   r="<<r_sb<<"+/-"<<er_sb<<endl;
	cout<<"L_b = "<<L_b<< ",   r="<<r_b<<"+/-"<<er_b<<endl;
	cout<<"L_sb_float = "<<L_sb_float<< ",   r="<<r_sb_float<<"+/-"<<er_sb_float<<endl;
	cout<<"L_b_float = "<<L_b_float<<endl;


	cout<<"number of independent sys sources = "<<cms->Get_max_uncorrelation()<<endl;
	return 1;
}
