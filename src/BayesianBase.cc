#include <math.h>
#include <float.h>
#include "time.h"
#include <cmath>
#include <stdio.h>
#include <algorithm>   // min(x,y)
#include <iostream>

// my include
#include "BayesianBase.h"
#include "Utilities.h"

//----------------need implementing a super good technique to do integration


using namespace std;
namespace lands{
	BayesianBase::BayesianBase(){
		_cms=0;
		fAlpha=0.05;
		fPrecision=1.e-3;
		fConfidenceLevel=1-fAlpha;
		_debug=0;
		_nexps_to_averageout_sys = 20000;
		_d=0;
		_xgl=0; _lwgl=0;
		_ngl=0;
		_norm=0; _stot=0; _btot=0;
	}
	BayesianBase::BayesianBase(double alpha, double precision){
		_cms=0;
		fAlpha=alpha;
		fPrecision=precision;
		fConfidenceLevel=1-fAlpha;
		_debug=0;
		_nexps_to_averageout_sys = 20000;
		_d=0;
		_xgl=0; _lwgl=0;
		_ngl=0;
		_norm=0; _stot=0; _btot=0;
	}
	BayesianBase::BayesianBase(CountingModel* cms, double alpha, double precision){
		_cms=cms;
		fAlpha=alpha;
		fPrecision=precision;
		fConfidenceLevel=1-fAlpha;
		_debug=0;
		_nexps_to_averageout_sys = 2000;
		_d=0;
		_nchannels = _cms->NumOfChannels();
		_xgl=0; _lwgl=0;
		_ngl=0;
		_norm=0; _stot=0; _btot=0;
	}
	void BayesianBase::SetModel(CountingModel *cms){
		//if(cms==_cms){cout<<" You are setting the same model ... "<<endl; return; }
		_cms=cms;
		_nchannels = _cms->NumOfChannels();
		_ngl=0;
	}
	BayesianBase::~BayesianBase(){
		_cms=0;
		for(int i=0; i<_vs.size(); i++){
			if(_vs[i]) delete [] _vs[i];
			if(_vb[i]) delete [] _vb[i];
		}
		_vs.clear();
		_vb.clear();
		if(_d) delete [] _d; 
		if(_xgl) delete [] _xgl;
		if(_lwgl) delete [] _lwgl;
		if(_stot) delete [] _stot;
		if(_btot) delete [] _btot;
	}
	double BayesianBase::Limit(CountingModel *cms){
		_cms=cms;	
		return Limit();
	}
	double BayesianBase::Limit(double alpha){
		fAlpha = alpha;
		fConfidenceLevel = 1-fAlpha;

		double ret=0;
		if(!_cms) {cout<<" model doesn't exist, exit"<<endl; exit(0);}
		if(_debug>=100) {
			_cms->Print(); 
		}

		// keep in memory
		bool bsys = _cms->IsUsingSystematicsErrors();
		int ntoys=_nexps_to_averageout_sys;

		// quick calc with 100 toys 
		if(!bsys) _nexps_to_averageout_sys=1;	
		else _nexps_to_averageout_sys=100;	

		GenToys();
		_norm = AverageIntegral(0);
		if(_debug) cout<<" norm  = "<<_norm<<endl;

		// first try  ---to get p<fAlpha
		double rlo= 0;
		double rmid = 10;
		double rhi=20;
		double p = AverageIntegral(rmid)/_norm;
		if(_debug) cout<<"p [ "<<rmid<<", inf ]"<<p<<endl;
		while(p>fAlpha) {
			rlo=rmid;
			rmid*=2;
			p = AverageIntegral(rmid)/_norm;
			if(_debug) cout<<"p [ "<<rmid<<", inf ] = "<<p<<endl;
		}
		rhi=rmid;
		
		// to converge at fAlpha+/-fPrecision
		while( fabs(p-fAlpha)/fAlpha > fPrecision ) {
			if(p<fAlpha){
				rhi=rmid;
				rmid=rlo+0.5*(rhi-rlo);				
				p = AverageIntegral(rmid)/_norm;
				if(_debug) cout<<"p [ "<<rmid<<", inf ] = "<<p<< "\t rlo="<<rlo<<" rhi="<<rhi<<endl;
			}else{
				rlo=rmid;
				rmid=rhi-0.5*(rhi-rlo);
				p = AverageIntegral(rmid)/_norm;
				if(_debug) cout<<"p [ "<<rmid<<", inf ] = "<<p<< "\t rlo="<<rlo<<" rhi="<<rhi<<endl;
			}	
		}

		ret = rmid;
		if(!bsys){
			_limit=ret;
			return ret;
		}


		if(_debug)cout<<"\t limit starts at = "<<ret<<endl;

		// resume the memory
		_cms->SetUseSystematicErrors(bsys);
		_nexps_to_averageout_sys=ntoys;
		GenToys();

		// try to get p<fAlpha 
		_norm = AverageIntegral(0);
		if(_debug) cout<<"\t ** norm (w/ sys) = "<<_norm<<" **"<<endl;

/*
		rmid = 2*ret;	
		p = AverageIntegral(rmid)/_norm;
		if(_debug) cout<<"p [ "<<rmid<<", inf ]"<<p<<endl;
		while(p>fAlpha) {
			rlo=rmid;
			rmid*=2;
			p = AverageIntegral(rmid)/_norm;
			if(_debug) cout<<"p [ "<<rmid<<", inf ] = "<<p<<endl;
		}
		rhi=rmid;

		// to converge at fAlpha+/-fPrecision
		while( fabs(p-fAlpha)/fAlpha > fPrecision ) {
			if(p<fAlpha){
				rhi=rmid;
				rmid=rlo+0.5*(rhi-rlo);				
				p = AverageIntegral(rmid)/_norm;
				if(_debug) cout<<"p [ "<<rmid<<", inf ] = "<<p<< "\t rlo="<<rlo<<" rhi="<<rhi<<endl;
				}else{
				rlo=rmid;
				rmid=rhi-0.5*(rhi-rlo);
				p = AverageIntegral(rmid)/_norm;
				if(_debug) cout<<"p [ "<<rmid<<", inf ] = "<<p<< "\t rlo="<<rlo<<" rhi="<<rhi<<endl;
				}	
				}
				ret = rmid;
 */
		// Joel Heinrich's approach to converge	
		const double eps=1.0e-6;
		double norm = _norm;
		double limit = rmid;
		double dl=limit, rpdf=0;
		double lo=0, hi=1e200;
		double beta = fConfidenceLevel;

		while(fabs(dl)>1.0e-4*limit) {
			const double pbeta =
				1-AverageIntegral(limit)/norm;
			rpdf = 1/Likelihood(limit);
			if(_debug) cout<<"p [ "<<limit<<", inf ] = "<<1-pbeta<<"   L("<<limit<<")= "<<1/rpdf<<endl;
			if(pbeta>beta) {
				hi=limit*(1+eps);
			} else {
				lo=limit*(1-eps);
			}
			dl = (pbeta-beta)*rpdf;
			if(limit-dl>=lo && limit-dl<=hi) {
				limit -= dl;
			} else {
				dl = limit - 0.5*(hi+lo);
				limit = 0.5*(hi+lo);
			}
		}
		ret=limit;

		_limit=ret;
		return ret;

	}
	double BayesianBase::LimitMeanInMath(){
		/*
		   if(_debug) _cms->Print();
		   if(_cms->NumOfChannels()!=1) {
		   cout<<"LimitMeanInMath is for single channel, exit"<<endl;
		   exit(0);
		   }
		   bool bsys = _cms->IsUsingSystematicsErrors();
		   if(_cms->IsUsingSystematicsErrors()) {
		   cout<<"LimitMeanInMath will run without systematics errors"<<endl;
		   _cms->SetUseSystematicErrors(false);
		   }
		   double minpoisson = 1.e-5; //configurable
		   double r_avr = 0.;
		   double nbkg = _cms->GetExpectedNumber(0, 1);
		   int nmax=(int)(nbkg*2+1);
		   double delta = Poisson(nmax,nbkg) - minpoisson;
		   while(delta>0){
		   nmax=int(2*nmax);
		   delta = Poisson(nmax,nbkg) - minpoisson;
		   }

		   for(int n=0; n<nmax; n++){
		   _cms->AddObservedData(0, n);

		   r_avr+=Limit()*Poisson(n,nbkg);	
		   }

		   _cms->SetUseSystematicErrors(bsys);
		   return r_avr;
		 */
	}
	void BayesianBase::GenToys(){
		if(_debug >= 100) _cms->Print();
		clock_t start_time=clock(), cur_time=clock(); // timing
		for(int i=0; i<_vs.size(); i++){
			if(_vs[i]) delete [] _vs[i];
			if(_vb[i]) delete [] _vb[i];
		}
		_vs.clear();
		_vb.clear();
		if(_d) delete [] _d; 
		if(_stot) delete [] _stot;
		if(_btot) delete [] _btot;

		_d = new double[_nchannels];
		_stot=new double[_nexps_to_averageout_sys];
		_btot=new double[_nexps_to_averageout_sys];

		double dtot=0;
		_logscale=0;
		for(int ch=0; ch<_nchannels; ch++){	
			_d[ch]=_cms->Get_v_data()[ch];
			dtot+=_d[ch];
			_logscale-=lgamma(_d[ch]+1);
		}
		if(dtot!=_dtot) _ngl=0;
		_dtot=dtot;

		if(_debug){ start_time=cur_time; cur_time=clock(); cout << "\t TIME_in_lgamma " << (cur_time - start_time) << " microsec\n";}
		if(_ngl<=0) { 
			_ngl = 1 + (int)dtot/2;
			if(_xgl) delete [] _xgl;
			if(_lwgl) delete [] _lwgl;
			_xgl = new double[_ngl];
			_lwgl = new double[_ngl];
			gausslaguerre(_xgl, _lwgl, _ngl, 0.0);
			if(_debug){ start_time=cur_time; cur_time=clock(); cout << "\t TIME_in_gausslaguerre " << (cur_time - start_time) << " microsec\n";}
		}

		vector< vector<double> > vv;
		for(int i=0; i<_nexps_to_averageout_sys; i++){
			vv = _cms->FluctuatedNumbers();
			double *s = new double[_nchannels];
			double *b = new double[_nchannels];
			_stot[i]=0; _btot[i]=0;
			for(int ch=0; ch<_nchannels; ch++){	
				double totbkg = 0; 
				for(int isamp=1; isamp<vv[ch].size(); isamp++){
					totbkg+=vv[ch][isamp];
				}
				s[ch]=vv[ch][0];
				b[ch]=totbkg; //FIXME just sum up the bkgs
				if(_debug>=100 && i<100 ){
					cout<<"ch = "<<ch<<" s="<<s[ch]<<" b="<<b[ch]<<" d="<<_d[ch]<<endl;
				}
				_stot[i]+=s[ch];
				_btot[i]+=b[ch];
			}
			_vs.push_back(s);
			_vb.push_back(b);
		}
		if(_debug){ start_time=cur_time; cur_time=clock(); cout << "\t TIME_in_gentoys " << (cur_time - start_time) << " microsec\n";}
	}
	double BayesianBase::AverageIntegral(double rlow ){
		if (rlow<0) {cout<< " integral starting from <0, exit" <<endl; exit(0);}
		double ret=0;
		for(int i=0; i<_nexps_to_averageout_sys; i++){
			ret += glintegral(rlow, i);
		}
		ret/=(double)_nexps_to_averageout_sys;
		return ret;
	}
	double BayesianBase::glintegral(double rlow, int iexps ) {
		int i,k;
		double sum=0, rstot;
		if(_stot[iexps]<=0){ cout<<"total signal <= 0 ,exit "<<endl; exit(0); }
		rstot=1./_stot[iexps];

		for(k=0;k<_ngl;++k) {
			const double xr = _xgl[k]*rstot + rlow;
			double t = -rlow * _stot[iexps]  - _btot[iexps] + _logscale , v;
			for(i=0;i<_nchannels;++i)
				if(_d[i]>0)
					t += _d[i] * log( xr*_vs[iexps][i] + _vb[iexps][i] );
			sum += v = exp(_lwgl[k]+t);
			if(v<DBL_EPSILON*sum) break;
		}
		// PRIOR flat
	//	sum *= rstot;
		return sum;
	}
	double BayesianBase::Likelihood(double r){
		if (r<0) {cout<< " r <0, exit" <<endl; exit(0);}
		if (_norm<=0)  {cout<< " _norm <=0, exit" <<endl; exit(0);}
		double ret=0;
		int i, c;
		double t;
		for(i=0; i<_nexps_to_averageout_sys; i++){
			t= -_btot[i] - r*_stot[i]+ _logscale;
			for(c=0; c<_nchannels; c++)
				if(_d[c]>0) 
					t += _d[c] * log( _vb[i][c] + r*_vs[i][c] );
			//if(prior==flat)
			//ret += exp(t);
			ret += _stot[i] * exp(t);
		}
		return ret/(_norm*_nexps_to_averageout_sys);

	}
	void BayesianBase::PosteriorPdf(int bins, double rmin, double rmax){
		_vr.clear(); _vp.clear();
		if(rmax<=0) rmax=2*_limit;	
		if(rmax<=rmin){
			cout<<"Warning rmax < rmin "<<endl;
			return;
		}
		if(rmin<0) rmin=0;
		if(bins<=0) bins=1;
		double step = (rmax-rmin)/(double)bins;
		for(; rmin<=rmax; rmin+=step){
			_vr.push_back(rmin);	
			_vp.push_back(Likelihood(rmin));		
		}
	}
};

