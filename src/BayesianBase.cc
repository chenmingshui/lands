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

#include "RooDataSet.h"

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
		_prior=flat;
		_NormReduction = 0;
		_preToys = 100;
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
		_prior=flat;
		_NormReduction = 0;
		_preToys = 100;
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
		_prior=flat;
		_NormReduction = 0;
		_preToys = 100;
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
		else _nexps_to_averageout_sys=_preToys;	

		GenToys();

		_NormReduction = EvaluateNormReduction();
		if(_debug) cout<<"  _NormReduction = "<<_NormReduction<<endl;

		_norm = AverageIntegral(0);
		if(_debug) cout<<(bsys?"w/ sys, ":"w/o sys,")<<" norm  = "<<_norm<<" , with my own approach to converge"<<endl;
		if(_norm ==0 || _norm>8.2e307) {
			cout<<"ERROR: the normalization of likelihood is very unlikely  = "<<(_norm==0?0:_norm)<<endl;
			cout<<"ERROR: we set the limit to 0"<<endl;
			_limit = 0;
			return 0;
		}

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

		
		if(_nexps_to_averageout_sys<=_preToys) { _limit=ret; return ret;}

		GenToys();

		// try to get p<fAlpha 
		_norm = AverageIntegral(0);
		if(_debug) cout<<"\t ** norm "<<(bsys?"w/ sys = ":"w/o sys =")<<_norm<<" ** with Joel Hinrich's approach"<<endl;
		if(_norm ==0 || _norm>8.2e307) {
			cout<<"WARNING: the normalization of likelihood is very unlikely  = "<<(_norm==0?0:_norm)<<endl;
			cout<<"WARNING: we set the limit to pre-calculated "<<ret<<endl;
			_limit = ret;
			return ret;
		}

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

		if(_debug>=100)cout<<"DELETEME 1"<<endl;

		_vs.clear();
		_vb.clear();
		if(_d) delete [] _d; 
		if(_stot) delete [] _stot;
		if(_btot) delete [] _btot;

		if(_debug>=100)cout<<"DELETEME 2"<<endl;
		_vNorms_forShapeChannels.clear();
		_vParams_forShapeChannels.clear();

		_d = new double[_nchannels];
		_stot=new double[_nexps_to_averageout_sys];
		_btot=new double[_nexps_to_averageout_sys];

		double dtot=0;
		_logscale=0;
		if(_debug>=100)cout<<"DELETEME 3"<<endl;
		for(int ch=0; ch<_nchannels; ch++){	
			_d[ch]=_cms->Get_v_data()[ch];
			dtot+=_d[ch];
			_logscale-=lgamma(_d[ch]+1);
		}
		if(_debug>=100)cout<<"DELETEME 4"<<endl;
		for(int ch=0; ch<_cms->Get_vv_pdfs().size(); ch++){
			if(_debug>=100)cout<<"v_pdfs_roodataset.size = "<<_cms->Get_v_pdfs_roodataset().size()<<endl;
			if(_cms->Get_v_pdfs_roodataset()[0]) _cms->Get_v_pdfs_roodataset()[0]->Print();
			else cout<<" 0th: not exist "<<endl;
			dtot+=(_cms->Get_v_pdfs_roodataset()[ch])->sumEntries();
			if(_debug>=100)cout<<"DELETEME 42"<<endl;
		}
		if(_debug>=100)cout<<"DELETEME 5"<<endl;
		if(dtot!=_dtot) _ngl=0;
		_dtot=dtot;

		if(_debug){ start_time=cur_time; cur_time=clock(); cout << "\t TIME_in_lgamma " << (cur_time - start_time) << " microsec\n";}
		if(_ngl<=0) { 
			if(_prior == flat || _prior==corr)_ngl = 1 + (int)dtot/2;
			//if(_prior == prior_1overSqrtS)_ngl = 1 + (int)((dtot-0.5)/2);
			if(_prior == prior_1overSqrtS)_ngl = 1 + (int)(dtot/2) + 1000; // +100  is for more accurate, otherwise for lower fluctuation, the results are not good 
			if(_xgl) delete [] _xgl;
			if(_lwgl) delete [] _lwgl;
			_xgl = new double[_ngl];
			_lwgl = new double[_ngl];
			gausslaguerre(_xgl, _lwgl, _ngl, 0.0);
			if(_debug){ start_time=cur_time; cur_time=clock(); cout << "\t TIME_in_gausslaguerre " << (cur_time - start_time) << " microsec\n";}
		}

		if(_debug>=10)  cout<<"_nexps_to_averageout_sys="<<_nexps_to_averageout_sys<<endl;
		vector< vector<double> > vv;
		vector<int> v_sigproc = _cms->Get_v_sigproc();
		for(int i=0; i<_nexps_to_averageout_sys; i++){
			if(_debug>=100)cout<<" _nexps_to_averageout_sys: "<<i<<endl;
			vv = _cms->FluctuatedNumbers();
			if(_debug>=100)cout<<" after fluctuation: "<<i<<endl;
			double *s = new double[_nchannels];
			double *b = new double[_nchannels];
			_stot[i]=0; _btot[i]=0;
			for(int ch=0; ch<_nchannels; ch++){	
				double totbkg = 0; 
				double totsig = 0; 
				for(int isamp=v_sigproc[ch]; isamp<vv[ch].size(); isamp++){
					totbkg+=vv[ch][isamp];
				}
				for(int isamp=0; isamp<v_sigproc[ch]; isamp++){
					totsig+=vv[ch][isamp];
				}
				s[ch]=totsig;
				b[ch]=totbkg; //FIXME just sum up the bkgs
				if(_debug>=100 && i<100 ){
					cout<<"ch = "<<ch<<" s="<<s[ch]<<" b="<<b[ch]<<" d="<<_d[ch]<<endl;
				}
				_stot[i]+=s[ch];
				_btot[i]+=b[ch];
			}
			_vNorms_forShapeChannels.push_back(_cms->Get_vv_pdfs_norm_varied());
			_vParams_forShapeChannels.push_back(_cms->Get_v_pdfs_floatParamsVaried());
			for(int ch=0; ch<_cms->Get_vv_pdfs().size(); ch++){
				for(int isamp=0; isamp<_cms->Get_vv_pdfs()[ch].size(); isamp++){
					if(isamp<(_cms->Get_v_pdfs_sigproc())[ch])
						_stot[i]+=(_cms->Get_vv_pdfs_norm_varied())[ch][isamp];
					else _btot[i]+= (_cms->Get_vv_pdfs_norm_varied())[ch][isamp];
				}
			}
			if(_debug>=100 && i<100 ){
				cout<<" toy  "<<i<<" stot="<<_stot[i]<<" btot="<<_btot[i]<<endl;
			}
			_vs.push_back(s);
			_vb.push_back(b);
		}
		if(_debug){ start_time=cur_time; cur_time=clock(); cout << "\t TIME_in_gentoys " << _nexps_to_averageout_sys <<" toys, "<< (cur_time - start_time) << " microsec\n";}
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


		VChannelVSample vvs, vvb; vvs.clear(); vvb.clear();
		for(k=0;k<_ngl;++k) {
			const double xr = _xgl[k]*rstot + rlow;
			double t = -rlow * _stot[iexps]  - _btot[iexps] + _logscale , v;
			for(i=0;i<_nchannels;++i)
				if(_d[i]>0)
					t += _d[i] * log( xr*_vs[iexps][i] + _vb[iexps][i] );
			t+=_cms->EvaluateGL(_vNorms_forShapeChannels[iexps], _vParams_forShapeChannels[iexps], xr, vvs, vvb);
			if(_prior == prior_1overSqrtS)t -= 0.5*log(_xgl[k] + rlow * _stot[iexps] );

			//cout<< " _lwgl["<<k<<"]="<< _lwgl[k] <<" + t="<<t <<" = "<<_lwgl[k]+t<<endl;
			sum += v = exp(_lwgl[k]+t - _NormReduction);
			//cout<< " exp of above = "<<v<<endl;
			if(v<DBL_EPSILON*sum) break;
		}
		if(_prior==flat )// PRIOR flat
			sum *= rstot;
		if(_prior==prior_1overSqrtS)// PRIOR flat
			sum *= sqrt(rstot);
		return sum;
	}
	double BayesianBase::Likelihood(double r){
		if (r<0) {cout<< " r <0, exit" <<endl; exit(0);}
		if (_norm<=0)  {cout<< " _norm <=0, i.e. likelihood norm of posterior pdf is very unlikeli. exit  " << _norm <<endl; exit(0);}
		double ret=0;
		int i, c;
		double t;
		for(i=0; i<_nexps_to_averageout_sys; i++){
			t= -_btot[i] - r*_stot[i]+ _logscale;
			for(c=0; c<_nchannels; c++)
				if(_d[c]>0) 
					t += _d[c] * log( _vb[i][c] + r*_vs[i][c] );

			VChannelVSample vvs, vvb; vvs.clear(); vvb.clear();
			t+=_cms->EvaluateGL(_vNorms_forShapeChannels[i], _vParams_forShapeChannels[i], r, vvs, vvb);
			//for(c=0; c<_cms->Get_vv_pdfs().size(); c++){
			//	t+=_cms->EvaluateGL(c, r);
			//}

			if(_prior==flat)
				ret += exp(t - _NormReduction);
			if(_prior==corr){
				ret += _stot[i] * exp(t);
				//ret +=  exp(t)/_stot[i];
			}
			if(_prior==prior_1overSqrtS){
				if(r>0) ret += exp(t)/sqrt(r);
			}
		}
		return ret/(_norm*_nexps_to_averageout_sys);

	}
	double BayesianBase::PosteriorPdf(int bins, double rmin, double rmax){
		_vr.clear(); _vp.clear();
		if(rmax<=0) rmax=2*_limit;	
		if(rmax<=rmin){
			cout<<"Warning rmax <= rmin "<<endl;
			return 0;
		}
		if(rmin<0) rmin=0;
		if(bins<=0) bins=1;
		double step = (rmax-rmin)/(double)bins;
		for(; rmin<=rmax; rmin+=step){
			_vr.push_back(rmin);	
			_vp.push_back(Likelihood(rmin));		
		}

		double delta_r = 0;
		// will return computational error on the r due to iteration cutout, the fPrecision,  in case of no systematics
		// from the posterior pdf curve, you can get  
		// delta_alpha = Likelihood(rlimit) * delta_r
		delta_r = 1./Likelihood(_limit) * fAlpha*fPrecision;
		return delta_r;
	}
	double BayesianBase::ErrorOnR_DueToPrecision(){
		// will return computational error on the r due to iteration cutout, the fPrecision,  in case of no systematics
		// from the posterior pdf curve, you can get  
		// delta_alpha = Likelihood(rlimit) * delta_r
		return  1./Likelihood(_limit) * fAlpha*fPrecision;
	}
	double BayesianBase::ErrorOnR_DueToFiniteToys(int ntrials, double& mean){
		double* rtmp = new double[ntrials];
		double tot = 0;
		for(int i=0; i<ntrials; i++){
			rtmp[i]=Limit();
			tot+=rtmp[i];
		}
		double meantmp = tot/(double)ntrials;
		double v = 0;
		for(int i=0; i<ntrials; i++){
			v+=( (rtmp[i]-meantmp)*(rtmp[i]-meantmp) );
		}
		delete [] rtmp;
		mean = meantmp;
		return sqrt(v/(double)ntrials);
	}
	int BayesianBase::EvaluateNormReduction(){
		double ret = 0 ;

		int iexps=0; double rlow = 0;
		int i,k;
		double sum=0, rstot;
		if(_stot[iexps]<=0){ cout<<"total signal <= 0 ,exit "<<endl; exit(0); }
		rstot=1./_stot[iexps];

		VChannelVSample vvs, vvb; vvs.clear(); vvb.clear();
		for(k=0;k<_ngl;++k) {
			const double xr = _xgl[k]*rstot + rlow;
			double t = -rlow * _stot[iexps]  - _btot[iexps] + _logscale , v;
			for(i=0;i<_nchannels;++i)
				if(_d[i]>0)
					t += _d[i] * log( xr*_vs[iexps][i] + _vb[iexps][i] );
			t+=_cms->EvaluateGL(_vNorms_forShapeChannels[iexps], _vParams_forShapeChannels[iexps], xr, vvs, vvb);
			if(_prior == prior_1overSqrtS)t -= 0.5*log(_xgl[k] + rlow * _stot[iexps] );
			
			double tmp =  _lwgl[k]+t; 
			if(tmp>ret) ret = tmp;

			tmp -= 500*(int)(tmp/500);
			sum += v = exp( tmp );
			if(v<DBL_EPSILON*sum) break;
		}
		
		ret = 500*(int)(ret/500);
		return int(ret);
	}
};


