#include <math.h>
#include <cmath>
#include <stdio.h>
#include <algorithm>   // min(x,y)
#include <iostream>

// my include
#include "BayesianLimitBase.h"
#include "Utilities.h"

//----------------need implementing a super good technique to do integration


using namespace std;
namespace lands{
	double R_CL(CountingModel *cms, double cl_alpha, double precision, double MinLikelihood, double rUpperBound, bool debug){
		int nbin = cms->NumOfChannels();
		double *par=new double[3*nbin];
		for(int i=0; i<nbin; i++){	
			par[3*i]=cms->Get_vv_exp_sigbkgs()[i][0];
			par[3*i+1]=cms->Get_vv_exp_sigbkgs()[i][1];
			par[3*i+2]=cms->Get_v_data()[i];
		}
		double ret = R_CL(par, cl_alpha, nbin, precision, MinLikelihood, rUpperBound, debug);
		delete [] par;
		return ret;
	}
	double R_CL(vector<double> vs, vector<double> vb, vector<double> vd, double cl_alpha, double precision, double MinLikelihood, double rUpperBound, bool debug){
		int nbin = (int)vs.size();
		double *par=new double[3*nbin];
		for(int i=0; i<nbin; i++){	
			par[3*i]=vs[i];
			par[3*i+1]=vb[i];
			par[3*i+2]=vd[i];
		}
		double ret=R_CL(par, cl_alpha, nbin, precision, MinLikelihood, rUpperBound, debug);
		delete [] par;
		return ret;
	}
	double R_CL(vector<double> vs, vector<double> vb, double cl_alpha, double precision, double MinLikelihood, double rUpperBound, bool debug){
		int nbin = (int)vs.size();
		double *par=new double[3*nbin];
		for(int i=0; i<nbin; i++){	
			par[3*i]=vs[i];
			par[3*i+1]=vb[i];
			par[3*i+2]=vb[i];
		}
		double ret=R_CL(par, cl_alpha, nbin, precision, MinLikelihood, rUpperBound, debug);
		delete [] par;
		return ret;
	}
	double R_CL(double nsig, double nbkg, double cl_alpha, double precision, double MinLikelihood, double rUpperBound, bool debug){
		double par[3]={nsig, nbkg, nbkg};
		return R_CL(par, cl_alpha, 1, precision, MinLikelihood, rUpperBound, debug);
	}
	double R_CL(double nsig, double nbkg, int ndat, double cl_alpha, double precision, double MinLikelihood, double rUpperBound, bool debug){
		double par[3]={nsig, nbkg, ndat};
		return R_CL(par, cl_alpha, 1, precision, MinLikelihood, rUpperBound, debug);
	}
	double R_CL(double *par, double cl_alpha, int nbins, double precision, double MinLikelihood, double rUpperBound, bool debug){
		if(debug){
			cout<<endl<<"************ starting debug Bayesian::R_CL ***********"<<endl;
			cout<<"Bayesian::R_CL nchannels="<<nbins<<"  s b d=";
			for(int i=0; i<nbins; i++){
				cout<<par[3*i]<<" "<<par[3*i+1]<<" "<<par[3*i+2]<<", ";
			}
			cout<<endl;
			cout<<"precision="<<precision<<" minimum likelihood="<<MinLikelihood<<" rUpperBound="<<rUpperBound<<endl;
		}
		double delta; 
		int npar = 3*nbins; //nbins or nchannnels
		double rmax = rUpperBound;  // r = sigma / sigma_SM

		double norm = myIntegral(P_n0_Given_brs, 0, rUpperBound, par, npar, 10e4);
		/*
		   if(steps = 10e5), it takes much longer time than 10e3
		 */
		//cout<<"count mememememeememememememem"<<endl;
		if(debug) cout<<"myIntegral norm="<<norm<<endl;
		if(debug) {
			double	norm2 = Integral(P_n0_Given_brs, 0, rUpperBound, par, npar, precision);
			cout<<"Integral norm P_n0_Given_brs="<<norm2<<endl;

			double step=0.01;
			cout<<"myIntegral---------: step="<<step<<endl;
			double tmp=0;
			for(double r=0; r<rUpperBound; r+=step){
				double rtmp[1]={r+step/2.};
				tmp+=P_n0_Given_brs(rtmp,par,npar)*step;	
			}
			cout<<"--------myIntegral="<<tmp<<endl;
			//cout<<"Integral(P_n0_Given_brs(double *r, double *par),0, rUpperBound, par)="<<Integral(P_n0_Given_brs, 0, rUpperBound, par)<<endl;
			if(fabs((norm2-tmp)/tmp) > 10e-2) 
				cout<<"NOTE******use myIntegral="<<tmp<<" not "<<norm2<<endl;
			norm2=tmp;
		}

		double x[1];
		x[0]=rmax;
		delta = (Likelihood_r(x,par,npar)/norm-MinLikelihood)/MinLikelihood;

		while(delta>0){
			rmax=2*rmax;
			x[0]=rmax;
			delta = (Likelihood_r(x,par,npar)/norm-MinLikelihood)/MinLikelihood;
		}
		if(debug) cout<<"trying to reach delta<=0: rmax="<<rmax<<" delta="<<delta<<endl;

		double x1=0, x2=rmax;
		//while (fabs(delta) > precision) {
		while (fabs(delta) > 1.e-4) {
			if(delta>0) x1=x[0];
			else x2=x[0];
			x[0]=(x1+x2)/2.;
			delta = (Likelihood_r(x,par,npar)/norm-MinLikelihood)/MinLikelihood;	  
		}
		if(debug) cout<<"trying to reach precision: x2="<<x2<<" delta="<<delta<<endl;
		if(debug) cout<<"Likelihood_r("<<x[0]<<",par,npar)="<<Likelihood_r(x,par,npar)<<endl;
		rmax = x[0];

		int i=0;
		if(rmax>1){
			while(rmax>10.){
				rmax/=10.;
				i++;
			}
			int ii=i;
			rmax = (int)(rmax+1)*::pow(10,ii);
		}else{
			while(rmax<1.){
				rmax*=10.;
				i++;
			}
			int ii=i;
			rmax=(int)(rmax+1)*::pow(10,-ii);
		}	

		if(debug) cout<<"rmax="<<rmax<<endl;

		//	double norm = Integral(Likelihood_r, 0.,rmax,par,npar);
		//	double norm =1;
		//	cout<<"norm="<<norm<<endl;
		double r1 = 0;
		double r2 = rmax;
		rmax = (r1 + r2)/2.;

		if(debug) cout<<"rmax="<<rmax<<endl;

		delta = myIntegral(Likelihood_r, 0,rmax,par,npar)/norm - cl_alpha;
		if(debug)cout<<"Integral(Likelihood_r,"<<0<<","<<rmax<<",par,"<<npar<<","<<precision<<")/"<<norm<<"="<<Integral(Likelihood_r, 0,rmax,par,npar,precision)/norm<<endl;
		if(debug) cout<<"delta=Integral-cl_alpha="<<delta<<endl;
		int numofintg=2;
		//while (fabs(delta) > precision) 
		while (fabs(delta) > 1.e-4) 
		{
			numofintg++;
			if (delta < 0) {
				r1 = rmax;
				rmax = (r1 + r2)/2.;
				delta += myIntegral(Likelihood_r, r1,rmax,par,npar)/norm;
				if(debug)cout<<"Integral(Likelihood_r,"<<r1<<","<<rmax<<",par,"<<npar<<","<<precision<<")="<<Integral(Likelihood_r, r1,rmax,par,npar,precision)/norm<<endl;
			}
			else{
				r2 = rmax;
				rmax = (r1 + r2)/2.;
				delta -= myIntegral(Likelihood_r, rmax, r2, par,npar)/norm;
				if(debug)cout<<"Integral(Likelihood_r,"<<rmax<<","<<r2<<",par,"<<npar<<","<<precision<<")="<<Integral(Likelihood_r, rmax,r2,par,npar,precision)/norm<<endl;
			}

			if(debug)	std::cout<<"rmax="<<rmax<<" delta="<<delta<<std::endl;
		} 
		if(debug)	cout<<"R_CL number_of_integration="<<numofintg<<endl;
		if(debug){
			cout<<"************ end debug Bayesian::R_CL ***********"<<endl;
		}
		return rmax;
	}

	double Likelihood_r(double *r, double *par, int npar){
		double result = 0;
		//------this is the super time consuming part if you calc norm factor everytime
		//------if it's same for several times, then get the norm factor outside and do normalization just by dividing a number.
		//------this was spotted and it saved me ass 
		//	double norm = Integral(P_n0_Given_brs, 0, rUpperBound, par, npar);
		//	cout<<"Likelihood_r_norm="<<norm<<endl;
		double norm =1;

		//	std::cout<<"norm="<<norm<<std::endl;
		result = P_n0_Given_brs(r,par, npar)/norm;
		//	result = P_n0_Given_brs(r,par, npar);
		return result;
	}
	double P_n0_Given_brs(double *r, double *par){
		return	Poisson(par[2],par[1]+r[0]*par[0]);
	}
	double P_n0_Given_brs(double *r, double *par,int npar){
		// par[0]=Sig; // expected signal within SM
		// par[1]=Bkg; // expected backgrounds
		// par[2]=N;   // observed totol number of events
		double result = 1;
		for(int i=0; i<(int)npar/3; i++){  // 4 channels
			// -- optimise it
			result *= Poisson(par[2+i*3],par[1+i*3]+r[0]*par[0+i*3]);
			//			cout<<"Poisson("<<par[2+i*3]<<","<<par[1+i*3]+r[0]*par[0+i*3]<<")="<<Poisson(par[2+i*3],par[1+i*3]+r[0]*par[0+i*3])<<endl;
		}
		return result;
	}
	double R_CL_avr(double nsig, double nbkg, double cl_alpha){
		double minpoisson = 1.e-5; //configurable
		double r_avr = 0.;
		int nmax=(int)(nbkg*2+1);
		double delta = Poisson(nmax,nbkg) - minpoisson;
		while(delta>0){
			nmax=int(2*nmax);
			delta = Poisson(nmax,nbkg) - minpoisson;
		}

		for(int n=0; n<nmax; n++){
			double par[3] = {nsig, nbkg, n};
			r_avr+=R_CL(par, cl_alpha, 1)*Poisson(n,nbkg);	
		}
		return r_avr;
	}
	double myIntegral(double (*fcn)(double *, double *, int), double a, double b, double *par, int npar, double steps){
		if(b==a) return 0;
		if(b<a)  {cerr<<"integral upbound < downbound; a="<<a<<" b="<<b<<endl; return 0;}
		double step=(b-a)/steps;
		double tmp=0;
		for(double r=a; r<b; r+=step){
			double rtmp[1]={r+step/2.};
			tmp+=P_n0_Given_brs(rtmp,par,npar)*step;	
		}
		return tmp;
	}
	};

