#include "BayesianLimitBase.h"
#include "Utilities.h"
#include "BayesianLimit.h"
#include <iostream>
#include <cstdio>
#include <math.h>
namespace lands{
BayesianLimit::BayesianLimit(CountingModel *cms, int nexps, CRandom *rnd){
	_nexps=nexps;
	_rdmGen=rnd;
	_meanLimit=0;
	_medianLimit=0;
	_1SigmaLow=0;
	_1SigmaHigh=0;
	_2SigmaLow=0;
	_2SigmaHigh=0;
	debug=0;

	_nchannels=cms->NumOfChannels();
	for(int n=0; n<_nchannels; n++){
		_nsig[n]=(cms->Get_vv_exp_sigbkgs())[n][0];
	}
	for(int n=0; n<_nchannels; n++){
		_nbkg[n]=(cms->Get_vv_exp_sigbkgs())[n][1]; // currently dealwith one background per channel  FIXME
	}
	_cms=cms;
}
BayesianLimit::~BayesianLimit()
{
	_rdmGen=NULL;
	v_rLimits.clear();
	v_rLimits_n.clear();
	v_rLimits_sorted.clear();
	v_cumulaP_sorted.clear();
	_cms=0;
}
void BayesianLimit::SetRdm(CRandom *rnd){
	_rdmGen=rnd;
}
void BayesianLimit::Run(){
	if(_nchannels==1) {
		SingleChannelRun(_nsig[0], _nbkg[0]);
		return;
	}

	v_rLimits.clear();
	v_rLimits_n.clear();
	v_rLimits_sorted.clear();
	v_cumulaP_sorted.clear();
	vector<double> v_rLimits_tmp;
	vector<int> v_rLimits_n_tmp;
	v_rLimits_tmp.clear();
	v_rLimits_n_tmp.clear();
	if(_nchannels>10000) {
		cout<<"Error:::::::Are you sure to run on bins/channels more then 10000 ?" << endl;
		cout<<"if not, please correct it; if so, please rebuild this package with some modification"<<endl;
		exit(0);
	}


	//	if(_RunOnExistEnsemble==false){
	_differentialR95s.clear();
	double *pa=new double[_nchannels*3];
	if(debug)cout<<"BayesianLimit::Run _exps="<<_nexps<<endl;
	for(int i=0; i<_nexps; i++){
		for(int nc=0; nc<_nchannels; nc++){
			// ----- *fLumi --> expected limits with 200/pb
			pa[0+3*nc]= _nsig[nc]; 
			pa[1+3*nc]= _nbkg[nc]; 
			pa[2+3*nc]= _rdmGen->Poisson(_nbkg[nc]); 
		}
		double limit=0;
		limit = R_CL(pa,0.95, _nchannels);      
		_differentialR95s.push_back(limit);
		if(debug) {
			cout<<endl<<"outcomes: "<<" ";
			for(int j=0; j<_nchannels; j++){
				cout<<pa[2+3*j];
			}
			cout<<"  it's limit="<<limit<<endl;
		}

		bool added = false;
		for(int j=0; j<(int)v_rLimits.size(); j++){
			if(limit==v_rLimits[j]) {
				v_rLimits_n[j]++;
				v_rLimits_n_tmp[j]++;
				added=true;
			}
		}
		if(added==false) {
			v_rLimits.push_back(limit);
			v_rLimits_n.push_back(1);
			v_rLimits_tmp.push_back(limit);
			v_rLimits_n_tmp.push_back(1);
		}
	}
	delete [] pa;
	SortAndCumulative(_differentialR95s, v_rLimits_sorted, v_cumulaP_sorted);
	GetBandsByLinearInterpolation(v_rLimits_sorted, v_cumulaP_sorted, _1SigmaLow, _1SigmaHigh, _2SigmaLow, _2SigmaHigh);

}

double BayesianLimit::Limit_mean(){
	// --  Mean Limit
	double sum = 0;
	for(int n=0; n<(int)v_cumulaP_sorted.size(); n++){
		if(n==0) sum+=v_rLimits_sorted[n]*v_cumulaP_sorted[n];
		else sum+=v_rLimits_sorted[n]*(v_cumulaP_sorted[n]-v_cumulaP_sorted[n-1]);
	}
	_meanLimit=sum;
	return _meanLimit;
}
double BayesianLimit::Limit_sigma(int nsigma){

	if(debug) {
		cout<<" In GetMedianLimit, print out the r and p"<<endl;
		for(int i=0; i<(int)v_rLimits_sorted.size(); i++)
			printf("r %5.2f p %3.2f\n",v_rLimits_sorted[i], v_cumulaP_sorted[i]);
	}
	if(nsigma==0) {
	return GetBandByFermiCurveInterpolation(v_rLimits_sorted, v_cumulaP_sorted, 0.5);
	}else if(nsigma==-2){ 
		return _2SigmaLow;
	}else if(nsigma==-1){ 
		return _1SigmaLow;
	}else if(nsigma==2){ 
		return _2SigmaHigh;
	}else if(nsigma==1){ 
		return _1SigmaHigh;
	}else{
		cout<<"Error: unknown sigma, exit(0)"<<endl; 
		exit(0);
	}
}
double BayesianLimit::LimitForExpectedBonly(){
	double *pa=new double[_nchannels*3];
	for(int nc=0; nc<_nchannels; nc++){
		pa[0+3*nc]= _nsig[nc]; 
		pa[1+3*nc]= _nbkg[nc]; 
		pa[2+3*nc]= _nbkg[nc]; 
	}
	double limit=0;
	limit = R_CL(pa,0.95, _nchannels);      
	return limit;
}
	vector<double> BayesianLimit::GetDifferentialLimits(){return _differentialR95s;}
	vector<double> BayesianLimit::GetvLimits(){return v_rLimits_sorted;}
	vector<double> BayesianLimit::GetvLimits_CP(){return v_cumulaP_sorted;}
void BayesianLimit::SingleChannelRun(double nsig, double nbkg){
	// shortcoming: if nbkg < 1, say nbkg = 0.02, then Poisson(0, 0.02) = 0.98... 
	//              no meaning for 1sigma, 2sigma
	/*
	   if(_nchannels!=1) {
	   cout<<"\n!!!!!!!SingleChannelRun is supposed to be run on single channel !!!!!!\n"<<endl;
	   return;
	   }
	 */
	v_rLimits_sorted.clear();
	v_cumulaP_sorted.clear();
	const double cumpMax = 0.999;
	const int nMax = 10000;
	int nActual = 0;
	double rn[nMax]; // r(n0|b,s)
	double pn[nMax]; // p(n<=n0|b)
	for(int n=0; n<nMax; n++){
		double par[3]={nsig,nbkg,n};
		rn[n]=R_CL(par, 0.95, 1);
		pn[n]=0;
		for(int j=0; j<=n; j++){
			pn[n]+=Poisson(j, nbkg);
		}
		//	cout<<"n="<<n<<" "<<rn[n]<<" "<< pn[n]<<endl;
		nActual++;

		v_rLimits_sorted.push_back(rn[n]);
		v_cumulaP_sorted.push_back(pn[n]);

		if(pn[n]>cumpMax) break;
	}

	GetBandsByFermiCurveInterpolation(v_rLimits_sorted, v_cumulaP_sorted, _1SigmaLow, _1SigmaHigh, _2SigmaLow, _2SigmaHigh);
}
};
