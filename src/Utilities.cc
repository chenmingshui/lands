/*
 * =====================================================================================
 *
 *       Filename:  Utilities.cc
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  09/19/2009 10:43:51 PM CEST
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Mingshui Chen (), Mingshui.Chen@cern.ch
 *        Company:  Univsity of Florida
 *
 * =====================================================================================
 */
#include "Utilities.h"
#include <math.h>
#include <float.h>
#include <iostream>
#include <cmath>
#include <limits> 

using std::cout;
using std::endl;
namespace lands{
	void SortAndCumulative(vector<double> q_b_tmp, vector<double> & qn, vector<double>& pn, bool down){
		int nexps=q_b_tmp.size(); double *x=new double[nexps]; for(int i=0;i<nexps;i++)x[i]=q_b_tmp[i];
		SortAndCumulative(x,nexps,qn,pn, down); delete []x;
	}
	void SortAndCumulative(double *q_b_tmp, int nexps, vector<double> & qn, vector<double>& pn, bool down){
		// ---- sort as Q increased
		qn.clear(); pn.clear();
		int *iq_b = new int[nexps];
		Sort(nexps, q_b_tmp, iq_b, down);  // increased Q  

		double q;
		for(int i=0; i<nexps; i++){
			q= q_b_tmp[iq_b[i]];
			if(i==0) {
				qn.push_back(q);
				pn.push_back(1);
			}
			else{
				if(q==qn[ (int)qn.size()-1 ]){
					pn[ (int)qn.size()-1 ]++;
				}
				else{
					qn.push_back(q);
					pn.push_back(1);
				}
			}
		}
		double totalp=0;
		for(int i=0; i<pn.size(); i++){
			totalp+=pn[i];
		}
	//	cout<<"\t SortAndCumulative: totalp_supposeE ="<<totalp<<endl;
		totalp=0;
		for(int i=0; i<(int)pn.size(); i++){
			pn[i]/=(double)nexps;
			totalp+=pn[i];
			if(i>0) pn[i]+=pn[i-1];
		}
	//	cout<<"\t SortAndCumulative: totalp_suppose1 ="<<totalp<<endl;
		totalp=0;
		for(int i=0; i<pn.size(); i++){
			totalp+=pn[i];
		}
	//	cout<<"\t SortAndCumulative: totalp_supposeH ="<<totalp<<endl;
		delete [] iq_b;
	}
	int GetBands(vector<double> & vx, double& _1SigmaLow, double& _1SigmaHigh, double& _2SigmaLow, double& _2SigmaHigh){
		_1SigmaLow=0; _1SigmaHigh=0; _2SigmaHigh=0; _2SigmaLow=0;

		int nexps = (int)vx.size();	
		double *dx = new double[nexps];
		for(int i=0; i<nexps; i++){
			dx[i]=vx[i];
			//		if(dx[i]==0) cout<<"GetBands dx[i]: "<<dx[i]<<endl;
		}
		GetBands(dx, nexps, _1SigmaLow, _1SigmaHigh, _2SigmaLow, _2SigmaHigh);
		//	for(int i=0; i<nexps; i++){
		//		if(dx[i]==0) cout<<"1 dx["<<i<<"]=0"<<endl;
		//	}
		delete [] dx;
		return 1;
	}
	int GetBands(double *dx, int nexps, double& _1SigmaLow, double& _1SigmaHigh, double& _2SigmaLow, double& _2SigmaHigh){
		_1SigmaLow=0; _1SigmaHigh=0; _2SigmaHigh=0; _2SigmaLow=0;
		vector<double> rn, pn;
		SortAndCumulative(dx, nexps, rn, pn);		
		/*
		   rn.clear();pn.clear();
		   int *ix = new int[nexps];
		//	for(int i=0; i<nexps; i++){
		//		if(dx[i]==0) cout<<"2 dx["<<i<<"]=0"<<endl;
		//	}
		Sort(nexps, dx, ix, 0);
		//	for(int i=0; i<nexps; i++){
		//		if(dx[i]==0) cout<<"3 dx["<<i<<"]=0"<<endl;
		//		if(ix[i]>=nexps) cout<<"ix["<<i<<"]="<<ix[i]<<"> nexps="<<nexps<<endl;
		//	}
		for(int i=0; i<nexps; i++){
		double limit = dx[ix[i]];
		//		if(limit==0) cout<<"GetBands limit: "<<limit<<endl;
		bool added = false;
		for(int j=0; j<(int)rn.size(); j++){
		if(limit==rn[j]) {
		pn[j]++;
		added=true;
		}
		}
		if(added==false) {
		rn.push_back(limit);
		pn.push_back(1);
		}
		}
		for(int i=0; i<(int)pn.size(); i++){
		pn[i]/=(double)nexps;
		if(i>0) pn[i]+=pn[i-1];
		}
		delete [] ix;
		 */
		return	GetBands(rn, pn,  _1SigmaLow, _1SigmaHigh, _2SigmaLow, _2SigmaHigh);
		/*	
			cout<<endl<<endl<<"==========================start debugging GetBands========================="<<endl;
			for(int i=0; i<(int)pn.size(); i++){
			cout<<"x="<<rn[i]<<"  p="<<pn[i]<<endl;
			}
			cout<<"==========================end debugging GetBands========================="<<endl;
		 */
		return 1;
	}
	int GetBands(vector<double> & rn, vector<double>& pn,  double& _1SigmaLow, double& _1SigmaHigh, double& _2SigmaLow, double& _2SigmaHigh){
		_1SigmaLow=-1; _1SigmaHigh=-1; _2SigmaHigh=-1; _2SigmaLow=-1;
		// =====rn is already sorted and pn is cummulative probability
		double _GreenBandLow = (1- 0.683)/2.; //1 sigma
		double _GreenBandHigh = 1 - (1- 0.683)/2.; //1 sigma
		double _YellowBandLow = (1- 0.955)/2.; //2 sigma
		double _YellowBandHigh = 1 - (1- 0.955)/2.; //2 sigma
		for(int i=0; i<(int)rn.size();i++){
			if(_1SigmaLow<0) {
				if(pn[i]>=_GreenBandLow){
					_1SigmaLow=rn[i];
				}
			}
			if(_2SigmaLow<0) {
				if(pn[i]>=_YellowBandLow){
					_2SigmaLow=rn[i];
				}
			}
			if(_1SigmaHigh<0) {
				if(pn[i]>=_GreenBandHigh){
					_1SigmaHigh=rn[i];
				}
			}
			if(_2SigmaHigh<0) {
				if(pn[i]>=_YellowBandHigh){
					_2SigmaHigh=rn[i];
				}
			}
		}
		int nActual = (int) rn.size();
		if((_1SigmaHigh<0)||(_1SigmaLow<0)||(_2SigmaHigh<0)||(_2SigmaLow<0)) cout<<"GetBands Error: somting wrong with pn[n-1]="<<pn[nActual-1]<<endl;
		/*
		// G-Green; Y-Yellow; H-High; L-Low;y
		int nGLL=-1, nGHL=-1, nYLL=-1, nYHL=-1;
		int nGLH=-1, nGHH=-1, nYLH=-1, nYHH=-1;
		int nActual = (int) rn.size();
		for(int n=0; n<nActual; n++){
		if(pn[n] <= _GreenBandLow) nGLL=n;
		if(pn[n] >= _GreenBandLow && nGLH==-1) nGLH=n; 

		if(pn[n] <= _YellowBandLow) nYLL=n;
		if(pn[n] >= _YellowBandLow && nYLH==-1) nYLH=n;

		if(pn[n] <= _GreenBandHigh) nGHL=n;
		if(pn[n] >= _GreenBandHigh && nGHH==-1 ) nGHH=n;

		if(pn[n] <= _YellowBandHigh) nYHL=n;
		if(pn[n] >= _YellowBandHigh && nYHH==-1 ) nYHH=n;
		}

		if(nGLL<0) _1SigmaLow=rn[0];
		else{
		_1SigmaLow=LinearInterpolation(rn[nGLL],pn[nGLL], rn[nGLH],pn[nGLH], _GreenBandLow);
		}

		if(nYLL<0) _2SigmaLow=rn[0];
		else{
		_2SigmaLow=LinearInterpolation(rn[nYLL],pn[nYLL],rn[nYLH],pn[nYLH], _YellowBandLow);
		}

		if(nGHL<0) _1SigmaHigh=rn[0];
		else{
		_1SigmaHigh=LinearInterpolation(rn[nGHL],pn[nGHL],rn[nGHH],pn[nGHH],_GreenBandHigh);
		}

		if(nYHL<0) _2SigmaHigh=rn[0];
		else{
		_2SigmaHigh=LinearInterpolation(rn[nYHL],pn[nYHL],rn[nYHH],pn[nYHH], _YellowBandHigh);
		}
		//	if(!_2SigmaLow || !_2SigmaHigh || !_1SigmaHigh || !_1SigmaLow)
		//	{
		//		cout<<"GetBands Error: -2s --> 2s: "<<_2SigmaLow<<" "<<_1SigmaLow<<" "<<_1SigmaHigh<<" "<<_2SigmaHigh<<endl;
		//		cout<<"GetBands rn[0]:"<<rn[0]<<endl;
		//	}
		 */
		return 1;
	}

	double GetBandByLinearInterpolation(vector<double>  rn, vector<double> pn,  double sigma){
		if(sigma>1) return -1;
		double result=-1;
		// =====rn is already sorted and pn is cummulative probability
		int nGLL=-1, nGLH=-1;
		int nActual = (int) rn.size();
		for(int n=0; n<nActual; n++){
			if(pn[n] <= sigma) nGLL=n;
			if(pn[n] >= sigma && nGLH==-1) nGLH=n; 
		}
		if(nGLL<0) result=rn[0];
		else{
			result=LinearInterpolation(rn[nGLL],pn[nGLL], rn[nGLH],pn[nGLH], sigma);
		}
		return result;

	}
	double GetBandByFermiCurveInterpolation(vector<double> rn, vector<double> pn,  double sigma){
		if(sigma>1) return -1;
		double result=-1;
		// =====rn is already sorted and pn is cummulative probability
		int nGLL=-1, nGLH=-1;
		int nActual = (int) rn.size();
		for(int n=0; n<nActual; n++){
			if(pn[n] <= sigma) nGLL=n;
			if(pn[n] >= sigma && nGLH==-1) nGLH=n; 
		}
		if(nGLL<0) result=rn[0];
		else{
			result=FCInterpolation(rn[nGLL],pn[nGLL], rn[nGLH],pn[nGLH], sigma);
		}
		return result;

	}
	int GetBandsByFeldmanCousins(vector<double>  rn, vector<double> pn,  double& _1SigmaLow, double& _1SigmaHigh, double& _2SigmaLow, double& _2SigmaHigh){
		_1SigmaLow=-1; _1SigmaHigh=-1; _2SigmaHigh=-1; _2SigmaLow=-1;
		// =====rn is already sorted and pn is cummulative probability
		// ---check if rn is sorted as increased
		// ---for the -2lnQ sortation, it should be decreased 


		int nActual = (int) rn.size();

		double onesigma = 0.683;
		double twosigma = 0.955;
		double distance=999999.;
		for(int n=0; n<nActual; n++){
			if( n==0 &&  pn[n]>=onesigma) {
				_1SigmaLow  = rn[n];  _1SigmaHigh = rn[n];  break;
			}

			for(int n2=n+1; n2<nActual; n2++){
				if(pn[n2]-(n>0?pn[n-1]:0) >= onesigma) {
					if(rn[n2] - rn[n] < distance )  {
						distance = rn[n2] - rn[n];
						_1SigmaLow = rn[n]; _1SigmaHigh = rn[n2];
					}
					break;
				}
			}

			if( pn[n]> 1-onesigma )  break;
		}

		distance=999999.;
		for(int n=0; n<nActual; n++){
			if( n==0 &&  pn[n]>=twosigma) {
				_2SigmaLow  = rn[n];  _2SigmaHigh = rn[n];  break;
			}

			for(int n2=n+1; n2<nActual; n2++){
				if(pn[n2]-(n>0?pn[n-1]:0) >= twosigma) {
					if(rn[n2] - rn[n] < distance )  {
						distance = rn[n2] - rn[n];
						_2SigmaLow = rn[n]; _2SigmaHigh = rn[n2];
					}
					break;
				}
			}

			if( pn[n]> 1-twosigma )  break;
		}


		return 1;
	}
	int GetBandsByFermiCurveInterpolation(vector<double>  rn, vector<double> pn,  double& _1SigmaLow, double& _1SigmaHigh, double& _2SigmaLow, double& _2SigmaHigh){
		_1SigmaLow=-1; _1SigmaHigh=-1; _2SigmaHigh=-1; _2SigmaLow=-1;
		// =====rn is already sorted and pn is cummulative probability
		// ---check if rn is sorted as increased
		// ---for the -2lnQ sortation, it should be decreased 
		//	for(int i=1; i<rn.size(); i++){
		//		if(rn[i]<rn[i-1]) cout<<"GetBandsByFermiCurveInterpolation --- input rn not sorted correctly: rn["<<i<<"]="<<rn[i]<<" < rn["<<i-1<<"]="<<rn[i-1]<<endl;
		//	}
		double _GreenBandLow = (1- 0.683)/2.; //1 sigma
		double _GreenBandHigh = 1 - (1- 0.683)/2.; //1 sigma
		double _YellowBandLow = (1- 0.955)/2.; //2 sigma
		double _YellowBandHigh = 1 - (1- 0.955)/2.; //2 sigma
		// G-Green; Y-Yellow; H-High; L-Low;y
		int nGLL=-1, nGHL=-1, nYLL=-1, nYHL=-1;
		int nGLH=-1, nGHH=-1, nYLH=-1, nYHH=-1;
		int nActual = (int) rn.size();
		for(int n=0; n<nActual; n++){
			if(pn[n] <= _GreenBandLow) nGLL=n;
			if(pn[n] >= _GreenBandLow && nGLH==-1) nGLH=n; 

			if(pn[n] <= _YellowBandLow) nYLL=n;
			if(pn[n] >= _YellowBandLow && nYLH==-1) nYLH=n;

			if(pn[n] <= _GreenBandHigh) nGHL=n;
			if(pn[n] >= _GreenBandHigh && nGHH==-1 ) nGHH=n;

			if(pn[n] <= _YellowBandHigh) nYHL=n;
			if(pn[n] >= _YellowBandHigh && nYHH==-1 ) nYHH=n;
		}

		if(nGLL<0) _1SigmaLow=rn[0];
		else{
			_1SigmaLow=FCInterpolation(rn[nGLL],pn[nGLL], rn[nGLH],pn[nGLH], _GreenBandLow);
		}

		if(nYLL<0) _2SigmaLow=rn[0];
		else{
			_2SigmaLow=FCInterpolation(rn[nYLL],pn[nYLL],rn[nYLH],pn[nYLH], _YellowBandLow);
		}

		if(nGHL<0) _1SigmaHigh=rn[0];
		else{
			_1SigmaHigh=FCInterpolation(rn[nGHL],pn[nGHL],rn[nGHH],pn[nGHH],_GreenBandHigh);
		}

		if(nYHL<0) _2SigmaHigh=rn[0];
		else{
			_2SigmaHigh=FCInterpolation(rn[nYHL],pn[nYHL],rn[nYHH],pn[nYHH], _YellowBandHigh);
		}
		//	if(!_2SigmaLow || !_2SigmaHigh || !_1SigmaHigh || !_1SigmaLow)
		//	{
		//		cout<<"GetBands Error: -2s --> 2s: "<<_2SigmaLow<<" "<<_1SigmaLow<<" "<<_1SigmaHigh<<" "<<_2SigmaHigh<<endl;
		//		cout<<"GetBands rn[0]:"<<rn[0]<<endl;
		//	}
		return 1;

	}
	int GetBandsByLinearInterpolation(vector<double>  rn, vector<double> pn,  double& _1SigmaLow, double& _1SigmaHigh, double& _2SigmaLow, double& _2SigmaHigh){
		_1SigmaLow=-1; _1SigmaHigh=-1; _2SigmaHigh=-1; _2SigmaLow=-1;
		// =====rn is already sorted and pn is cummulative probability
		double _GreenBandLow = (1- 0.683)/2.; //1 sigma
		double _GreenBandHigh = 1 - (1- 0.683)/2.; //1 sigma
		double _YellowBandLow = (1- 0.955)/2.; //2 sigma
		double _YellowBandHigh = 1 - (1- 0.955)/2.; //2 sigma
		// G-Green; Y-Yellow; H-High; L-Low;y
		int nGLL=-1, nGHL=-1, nYLL=-1, nYHL=-1;
		int nGLH=-1, nGHH=-1, nYLH=-1, nYHH=-1;
		int nActual = (int) rn.size();
		for(int n=0; n<nActual; n++){
			if(pn[n] <= _GreenBandLow) nGLL=n;
			if(pn[n] >= _GreenBandLow && nGLH==-1) nGLH=n; 

			if(pn[n] <= _YellowBandLow) nYLL=n;
			if(pn[n] >= _YellowBandLow && nYLH==-1) nYLH=n;

			if(pn[n] <= _GreenBandHigh) nGHL=n;
			if(pn[n] >= _GreenBandHigh && nGHH==-1 ) nGHH=n;

			if(pn[n] <= _YellowBandHigh) nYHL=n;
			if(pn[n] >= _YellowBandHigh && nYHH==-1 ) nYHH=n;
		}

		if(nGLL<0) _1SigmaLow=rn[0];
		else{
			_1SigmaLow=LinearInterpolation(rn[nGLL],pn[nGLL], rn[nGLH],pn[nGLH], _GreenBandLow);
		}

		if(nYLL<0) _2SigmaLow=rn[0];
		else{
			_2SigmaLow=LinearInterpolation(rn[nYLL],pn[nYLL],rn[nYLH],pn[nYLH], _YellowBandLow);
		}

		if(nGHL<0) _1SigmaHigh=rn[0];
		else{
			_1SigmaHigh=LinearInterpolation(rn[nGHL],pn[nGHL],rn[nGHH],pn[nGHH],_GreenBandHigh);
		}

		if(nYHL<0) _2SigmaHigh=rn[0];
		else{
			_2SigmaHigh=LinearInterpolation(rn[nYHL],pn[nYHL],rn[nYHH],pn[nYHH], _YellowBandHigh);
		}
		//	if(!_2SigmaLow || !_2SigmaHigh || !_1SigmaHigh || !_1SigmaLow)
		//	{
		//		cout<<"GetBands Error: -2s --> 2s: "<<_2SigmaLow<<" "<<_1SigmaLow<<" "<<_1SigmaHigh<<" "<<_2SigmaHigh<<endl;
		//		cout<<"GetBands rn[0]:"<<rn[0]<<endl;
		//	}
		return 1;

	}
	//double lgamma(double x); // decalare double instance here
	int GetBandsByNoInterpolation(vector<double>  rn, vector<double> pn,  double& _1SigmaLow, double& _1SigmaHigh,  double& _2SigmaLow, double& _2SigmaHigh, double& median){
		_1SigmaLow=0; _1SigmaHigh=0; _2SigmaHigh=0; _2SigmaLow=0; median=0;
		// =====rn is already sorted and pn is cummulative probability
		double _GreenBandLow = (1- 0.683)/2.; //1 sigma
		double _GreenBandHigh = 1 - (1- 0.683)/2.; //1 sigma
		double _YellowBandLow = (1- 0.955)/2.; //2 sigma
		double _YellowBandHigh = 1 - (1- 0.955)/2.; //2 sigma
		bool brmedian2=false, brm1s2=false, brm2s2=false, brp1s2=false, brp2s2=false;	
		for(int i=0; i<pn.size(); i++){
			if(pn[i]>=_GreenBandLow && !brm1s2) { _1SigmaLow= rn[i]; brm1s2 = true; } 
			if(pn[i]>=_GreenBandHigh && !brp1s2) { _1SigmaHigh= rn[i]; brp1s2 = true; } 
			if(pn[i]>=_YellowBandLow && !brm2s2) { _2SigmaLow= rn[i]; brm2s2 = true; } 
			if(pn[i]>=_YellowBandHigh && !brp2s2) { _2SigmaHigh= rn[i]; brp2s2 = true; } 
			if(pn[i]>=0.5 && !brmedian2) { median= rn[i]; brmedian2 = true; } 
		}
		return 1;
	}
	double IntegralSum(double (*fcn)(double *, double *, int), double a, double b, double *par, int npar){
		// very dirty 
		double result = 0;
		if(a==b) return 0;
		if(a>b) {cout<<"ERROR: integral region a > b, it should be a <= b"<<endl; return 0;} 

		int nbins = 1000; //configurable
		double dx = (b-a)/nbins;

		double xx[1];

		double f = 0;

		double sdx = 0; 
		for(int i=0; i<nbins; i++){
			xx[0]=a+dx/2.;
			f = (*fcn)(xx, par, npar);
			if(f<0) f=-f;
			sdx=f*dx;
			result+=sdx;
		}

		//cout<<"integral "<<a<<" to " <<b<< " = "<<result<<endl;
		return result;

	}
	double Integral(double (*fcn)(double *, double *), double a, double b, double *par, double epsilon){
		//	cout<<"UtilitiesIntegral1"<<endl;
		// - - - from ROOT

		double fEpsilon = epsilon; // configurable  

		bool  fgAbsValue = true;

		const double kHF = 0.5;
		const double kCST = 5./1000;

		double x[12] = { 0.96028985649753623,  0.79666647741362674,
			0.52553240991632899,  0.18343464249564980,
			0.98940093499164993,  0.94457502307323258,
			0.86563120238783174,  0.75540440835500303,
			0.61787624440264375,  0.45801677765722739,
			0.28160355077925891,  0.09501250983763744};

		double w[12] = { 0.10122853629037626,  0.22238103445337447,
			0.31370664587788729,  0.36268378337836198,
			0.02715245941175409,  0.06225352393864789,
			0.09515851168249278,  0.12462897125553387,
			0.14959598881657673,  0.16915651939500254,
			0.18260341504492359,  0.18945061045506850};

		double h, aconst, bb, aa, c1, c2, u, s8, s16, f1, f2;
		double xx[1];
		int i;

		/*   if ( fFunction == 0 )
		     {
		     MATH_ERROR_MSG("ROOT::Math::GausIntegratorOneDim", "A function must be set first!");
		     return 0.0;
		     }
		 */
		h = 0;
		if (b == a) return h;
		if(b>a)aconst = kCST/(b-a);
		else aconst = - kCST/(b-a);

		bb = a;
CASE1:
		aa = bb;
		bb = b;
CASE2:
		c1 = kHF*(bb+aa);
		c2 = kHF*(bb-aa);
		s8 = 0;
		for (i=0;i<4;i++) {
			u     = c2*x[i];
			xx[0] = c1+u;
			//      f1    = (*fFunction)(xx);
			// --Mingshui
			//cout<<xx[0]<<" "<<par[0]<<" "<<par[1]<<endl;
			f1 = (*fcn)(xx, par);


			if (fgAbsValue) {if(f1<0)f1 = -f1;}
			xx[0] = c1-u;
			//f2    = (*fFunction) (xx);
			// --Mingshui
			f2 = (*fcn)(xx, par);

			if (fgAbsValue) {if(f2<0)f2 = -f2;}
			s8   += w[i]*(f1 + f2);
		}
		s16 = 0;
		for (i=4;i<12;i++) {
			u     = c2*x[i];
			xx[0] = c1+u;
			//f1    = (*fFunction) (xx);
			// --Mingshui
			f1 = (*fcn)(xx, par);

			if (fgAbsValue) {if(f1<0)f1 = -f1;}
			xx[0] = c1-u;
			//f2    = (*fFunction) (xx);
			// --Mingshui
			f2 = (*fcn)(xx, par);

			if (fgAbsValue) {if(f2<0)f2 = -f2;}
			s16  += w[i]*(f1 + f2);
		}
		s16 = c2*s16;
		double s16_tmp = s16;
		double s16_c2s8_tmp = s16-c2*s8;
		double c2_tmp=c2;


		if(s16_c2s8_tmp<0) s16_c2s8_tmp= -s16_c2s8_tmp;
		if(s16_tmp<0) s16_tmp=-s16_tmp;
		if(c2_tmp<0) c2_tmp=-c2_tmp;

		if (s16_c2s8_tmp <= fEpsilon*(1. + s16_tmp)) {
			h += s16;
			if(bb != b) goto CASE1;
		} else {
			bb = c1;
			if(1. + aconst*c2_tmp != 1) goto CASE2;
			h = s8;  //this is a crude approximation (cernlib function returned 0 !)
		}

		//  fUsedOnce = true;
		//  fLastResult = h;
		//  fLastError = std::abs(s16-c2*s8);

		return h;

	}

	double FermiCurve(double *x, double *par){
		if(par[1]==0) return 0;
		double result = 1./(1+exp((par[0]-x[0])/par[1]));
		return result;
	}

	void erase(vector<double>& v, int pos){
		int last_pos = (int)v.size()-1;
		double tmp1=v[last_pos];
		double tmp2=v[pos];
		v[pos]+=tmp1;
		v[pos]-=tmp2;
		//v[pos]=v[last_pos];
		v.pop_back();
	}
	void erase(vector<int>& v, int pos){
		int last_pos = v.size()-1;
		int tmp1= v[last_pos];
		int tmp2= v[pos];
		v[pos]+=tmp1;
		v[pos]-=tmp2;
		//v[pos]=v[last_pos];
		v.pop_back();
	}
	double LogLinearInterpolation(double x1, double y1, double x2, double y2, double y){
		if(y1<=0 || y2<=0 || y<=0 ) return 0;
		double x=0;
		if(x1==x2 || y1==y2) x=x1;
		else {
			x= log(y2/y)*(x1-x2)/log(y2/y1) + x2;
		}	
		return x;
	}
	double LogLinearInterpolationErr(double x1, double y1, double ey1, double x2, double y2, double ey2, double y){
		double x =    LogLinearInterpolation(x1, y1, x2, y2, y);
		double xmax = LogLinearInterpolation(x1, y1+ey1, x2, y2+ey2, y);
		double xmin = LogLinearInterpolation(x1, y1>ey1?(y1-ey1): y1, x2, y2>ey2?(y2-ey2):y2, y);
		double ex =  fabs(xmax - x);
		if(fabs(xmax-x) < fabs(xmin-x)) ex = fabs(xmin - x);
		return ex;
	}
	double LinearInterpolation(double x1, double y1, double x2, double y2, double y){
		double x=0;
		if(x1==x2 || y1==y2) x=x1;
		else {
			//if(y1==0 || y2==0) return 0; //FIXME
			double a=(y2-y1)/(x2-x1);
			double b=y2-a*x2;
			x=(y-b)/a;
		}	
		return x;
	}
	double FCInterpolation(double x1, double y1, double x2, double y2, double y){
		// fermi function:  y(x) = 1/(1+exp(-((x-x0)/tau))) 
		// --- detail in Andrey's note
		double x=0;
		if(x1==x2 || y1==y2) x=x1;
		else {
			if(y1==0 || y2==0) return 0;
			double lambda = log((1-y)/y);
			double lambda1 = log((1-y1)/y1);
			double lambda2 = log((1-y2)/y2);
			double tau = -(x2-x1)/(lambda2-lambda1);
			double x0  = x1+lambda1*tau; 
			x = x0-lambda*tau;
		}	
		return x;
	}

	double Integral(double (*fcn)(double *, double *, int), double a, double b, double *par, int npar, double epsilon){

		// - - - from ROOT
		//   based on original CERNLIB routine DGAUSS by Sigfried Kolbig
		//   converted to C++ by Rene Brun

		double fEpsilon = epsilon; // configurable  

		bool  fgAbsValue = true;

		const double kHF = 0.5;
		const double kCST = 5./1000;

		double x[12] = { 0.96028985649753623,  0.79666647741362674,
			0.52553240991632899,  0.18343464249564980,
			0.98940093499164993,  0.94457502307323258,
			0.86563120238783174,  0.75540440835500303,
			0.61787624440264375,  0.45801677765722739,
			0.28160355077925891,  0.09501250983763744};

		double w[12] = { 0.10122853629037626,  0.22238103445337447,
			0.31370664587788729,  0.36268378337836198,
			0.02715245941175409,  0.06225352393864789,
			0.09515851168249278,  0.12462897125553387,
			0.14959598881657673,  0.16915651939500254,
			0.18260341504492359,  0.18945061045506850};

		double h, aconst, bb, aa, c1, c2, u, s8, s16, f1, f2;
		double xx[1];
		int i;

		/*   if ( fFunction == 0 )
		     {
		     MATH_ERROR_MSG("ROOT::Math::GausIntegratorOneDim", "A function must be set first!");
		     return 0.0;
		     }
		 */
		h = 0;
		if (b == a) return h;
		if(b>a)aconst = kCST/(b-a);
		else aconst = - kCST/(b-a);

		bb = a;
CASE3:
		aa = bb;
		bb = b;
CASE4:
		c1 = kHF*(bb+aa);
		c2 = kHF*(bb-aa);
		s8 = 0;
		for (i=0;i<4;i++) {
			u     = c2*x[i];
			xx[0] = c1+u;
			//      f1    = (*fFunction)(xx);
			// --Mingshui
			f1 = (*fcn)(xx, par, npar);

			if (fgAbsValue) {if(f1<0)f1 = -f1;}
			xx[0] = c1-u;
			//f2    = (*fFunction) (xx);
			// --Mingshui
			f2 = (*fcn)(xx, par, npar);

			if (fgAbsValue) {if(f2<0)f2 = -f2;}
			s8   += w[i]*(f1 + f2);
		}
		s16 = 0;
		for (i=4;i<12;i++) {
			u     = c2*x[i];
			xx[0] = c1+u;
			//f1    = (*fFunction) (xx);
			// --Mingshui
			f1 = (*fcn)(xx, par,npar);

			if (fgAbsValue) {if(f1<0)f1 = -f1;}
			xx[0] = c1-u;
			//f2    = (*fFunction) (xx);
			// --Mingshui
			f2 = (*fcn)(xx, par, npar);

			if (fgAbsValue) {if(f2<0)f2 = -f2;}
			s16  += w[i]*(f1 + f2);
		}
		s16 = c2*s16;
		double s16_tmp = s16;
		double s16_c2s8_tmp = s16-c2*s8;
		double c2_tmp=c2;


		if(s16_c2s8_tmp<0) s16_c2s8_tmp= -s16_c2s8_tmp;
		if(s16_tmp<0) s16_tmp=-s16_tmp;
		if(c2_tmp<0) c2_tmp=-c2_tmp;

		if (s16_c2s8_tmp <= fEpsilon*(1. + s16_tmp)) {
			h += s16;
			if(bb != b) goto CASE3;
		} else {
			bb = c1;
			if(1. + aconst*c2_tmp != 1) goto CASE4;
			h = s8;  //this is a crude approximation (cernlib function returned 0 !)
		}

		//  fUsedOnce = true;
		//  fLastResult = h;
		//  fLastError = std::abs(s16-c2*s8);

		//	cout<<"UtilitiesIntegral2"<<endl;
		return h;

	}

	double Poisson(double *x, double *par){
		return Poisson(x[0],par[0]);
	}

	double Poisson(double x, double par)
	{
		// - - - -From ROOT TMath,    exp, log, lgamma
		// compute the Poisson distribution function for (x,par)
		// The Poisson PDF is implemented by means of Euler's Gamma-function
		// (for the factorial), so for all integer arguments it is correct.
		// BUT for non-integer values it IS NOT equal to the Poisson distribution.
		// see TMath::PoissonI to get a non-smooth function.
		// Note that for large values of par, it is better to call
		//     TMath::Gaus(x,par,sqrt(par),kTRUE)

		if(par<=0) return 0;
		if (x<0)
			return 0;
		else if (x == 0.0)
			return 1./exp(par);
		else {
			double lnpoisson = x*log(par)-par-lgamma(x+1.);
			return exp(lnpoisson);
		}
		// An alternative strategy is to transition to a Gaussian approximation for
		// large values of par ...
		//   else {
		//     return Gaus(x,par,Sqrt(par),kTRUE);
		//   }
	}
	double GetMeanOfSortedXwithProb(vector<double> vx, vector<double> vp ){
		// vx is supposed to be sorted, if not, you are in trouble 
		// FIXME check here if it's sorted
		// and vx.size = vp.size
		// here p is cummulative ...
		double mean=0;
		if(vx.size()==vp.size()){
			for(int i=0; i<vx.size(); i++){
				if(i==0){
					mean+=vx[i]*vp[i];
				}
				else {
					mean+=vx[i]*(vp[i]-vp[i-1]);
				}
			}	
		}
		return mean; 
	}
	double TruncatedGaussianPdf(double *x, double *par){
		double pi_tmp=3.14159265358979323846;
		double ret;
		if(par[0]<=0)return 0; // par[0]=sigma, should be > 0
		if(x[0]<-1) return 0;	
		else{
			ret=2./( 1+erf( 1./sqrt(2.)/par[0] ) )/( sqrt(2*pi_tmp) ) * exp(-x[0]*x[0]/2./par[0]/par[0]);	
		}
		return ret;
	}


	long nto1(int i){
		if(i<=1) return 1;
		else{
			long ret=1;
			for(int j=i; j>0; j--){
				ret*=j;
			}
			return ret;
		} 
		return 1;
	}
	double CLs_Analytics(double s, double b){
		double d=b;
		return CLs_Analytics(s,b,d);
	}
	double CLs_Analytics(double s, double b, int d){
		//	+++++++++++++++++ CL
		double cls3 = 0;
		double tmpsb=0;
		double tmpb=0;
		for(int i=0; i<=d; i++){
			double p = (double)nto1(i);
			tmpsb+=exp(-s-b)*pow(b+s,i)/p;
			tmpb+=exp(-b)*pow(b,i)/p;
		}
		cls3=tmpsb/tmpb;
		return cls3;
	}
	double CLs_Analytics(double s, double b, double d){
		double cls3=0;
		double tmpsb=0;
		double tmpb=0;
		double parsb[1]; parsb[0]=s+b;
		double parb[1]; parb[0]=b;
		tmpsb = Integral(Poisson, 0, d, parsb );	
		tmpb = Integral(Poisson, 0, d, parb );
		cls3=tmpsb/tmpb;
		return cls3;
	}
	vector< pair<double, double> > m2lnQ_b_Analytics(double s, double b){
		vector< pair<double, double> > tmp; tmp.clear();
		if(b==0) {cout<<"Error:  b=0"<<endl; return tmp; }
		double logsbb=log(1+s/b);
		double m2lnQ=0;
		for(int n=0; n<10000; n++){
			m2lnQ=-2*(-s+n*logsbb);
			double p=Poisson(n,b);
			tmp.push_back(pair<double,double>(m2lnQ,p));
			if(n>b && p<1.e-6) break; 
		}
		return tmp;
	}
	vector< pair<double, double> > m2lnQ_sb_Analytics(double s, double b){
		vector< pair<double, double> > tmp; tmp.clear();
		if(b==0) {cout<<"Error:  b=0"<<endl; return tmp; }
		double logsbb=log(1+s/b);
		double m2lnQ=0;
		for(int n=0; n<10000; n++){
			m2lnQ=-2*(-s+n*logsbb);
			double p=Poisson(n,b+s);
			tmp.push_back(pair<double,double>(m2lnQ,p));
			if(n>(b+s) && p<1.e-6) break; 
		}
		return tmp;
	}
	double m2lnQ(double s,double b, double d){
		if(s<0 || b<=0 || d<0) {
			cout<<"Error::::::::::::"<<endl;
			cout<<"s="<<s<<" b="<<b<<" d="<<d<<endl;
			cout<<"please make sure s>=0; b>0; d>=0"<<endl;
			return 0;
		}
		return -2*(-s+d*log(1+s/b));
	}

	double fP_n0_Given_brs(double *r, double *par,int npar){
		double result = 1;
		for(int i=0; i<(int)npar/3; i++){  
			result *= Poisson(par[2+i*3],par[1+i*3]+r[0]*par[0+i*3]);
		}
		return result;
	}

	/*
	   Joel Heinrich
	   February 10 2005

	   Returns Gauss-Laguerre quadrature abscissas and log(weights) which can
	   be used to approximate

	   integral u=0 to infinity pow(u,alpha)*exp(-u)*f(u) du
	   as
	   sum k=0 to n-1  exp(lw[k])*f(x[k])

	   or equivalently

	   sum k=0 to n-1  exp(lw[k]+log(f(x[k])))

	   The quadrature is exact for polynomial f of degree 2n-1 or less.

	 */
	void gausslaguerre(double x[],double lw[],int n,double alpha){
		const int nshift = 20;
		const double shift = 1<<nshift, rshift=1/shift;
		int i;
		double z=0;

		for(i=0;i<n;++i) {
			//cout<<"DELETEME in gausslaguerre "<<i<<" of n="<<n<<endl;
			int j=0, k=2, nscale=0;
			double dz=0.0, p1=0, p2=0;
			if(i==0) {
				z=(1.0+alpha)*(3.0+0.92*alpha)/(1.0+2.4*n+1.8*alpha);
			} else if(i==1) {
				z += (15.0+6.25*alpha)/(1.0+2.5*n+0.9*alpha);
			} else if(i==2) {
				const double ai=i-1;
				z += ( (1.0+2.55*ai)/(1.9*ai) + 1.26*ai*alpha/(1.0+3.5*ai) )*
					(z-x[i-2])/(1.0+0.3*alpha);
			} else if(i==3) {
				z = 3.0*(x[2]-x[1])+x[0];
			} else if(i==4) {
				z = 4.0*x[3] - 6.0*x[2] + 4.0*x[1] - x[0];
			} else if(i==5) {
				z = 5.0*x[4] - 10.0*x[3] + 10.0*x[2] - 5.0*x[1] + x[0];
			} else {
				z = 6.0*x[i-1] - 15.0*x[i-2] + 20.0*x[i-3] -
					15.0*x[i-4] + 6.0*x[i-5] - x[i-6];
			}
			while(k>0) {
				// when n >=7175,  it goes into infinite loop
				p1=1;
				p2=0;
				nscale=0;
				z -= dz;
				for(j=1;j<=n;++j){
					const double p3=p2;
					p2=p1;
					p1=((2*j-1+alpha-z)*p2 - (j-1+alpha)*p3)/j;
					if(fabs(p2)>shift) {
						++nscale;
						p1 *= rshift;
						p2 *= rshift;
					}
				}
				dz = p1*z/(n*p1-(n+alpha)*p2);
				//cout<<"DELETEME dz="<<dz<<" z="<<z<<endl;
				if(fabs(dz)<1.0e-10*z)--k;
				//cout<< " DELETEME k "<<k<< "  dz="<<dz<<"  z="<<z<<endl;
			}
			x[i]=z;
			lw[i] = log(z/(p2*p2)) - 2*nshift*nscale*M_LN2 ;
		}

		{
			double t = 0.0;
			for(i=n-1;i>=0;--i)
				t += exp(lw[i]);
			t = lgamma(alpha+1)-log(t);
			for(i=0;i<n;++i)
				lw[i] += t;
		}

		return;
	}

	namespace Cephes{
		/*
		 * calculates a value of a polynomial of the form:
		 * a[0]x^N+a[1]x^(N-1) + ... + a[N]
		 */
		double Polynomialeval(double x, double* a, unsigned int N)
		{
			if (N==0) return a[0];
			else
			{
				double pom = a[0];
				for (unsigned int i=1; i <= N; i++)
					pom = pom *x + a[i];
				return pom;
			}
		}

		/*
		 * calculates a value of a polynomial of the form:
		 * x^N+a[0]x^(N-1) + ... + a[N-1]
		 */
		double Polynomial1eval(double x, double* a, unsigned int N)
		{
			if (N==0) return a[0];
			else
			{
				double pom = x + a[0];
				for (unsigned int i=1; i < N; i++)
					pom = pom *x + a[i];
				return pom;
			}
		}


		// inverse of gamma and beta from Cephes library
		// see:  http://www.netlib.org/cephes
		// 
		// Copyright 1985, 1987, 2000 by Stephen L. Moshier

		/*
		 *
		 *      Inverse of Normal distribution function
		 *
		 *
		 *
		 * SYNOPSIS:
		 *
		 * double x, y, ndtri();
		 *
		 * x = ndtri( y );
		 *
		 *
		 *
		 * DESCRIPTION:
		 *
		 * Returns the argument, x, for which the area under the
		 * Gaussian probability density function (integrated from
		 * minus infinity to x) is equal to y.
		 *
		 *
		 * For small arguments 0 < y < exp(-2), the program computes
		 * z = sqrt( -2.0 * log(y) );  then the approximation is
		 * x = z - log(z)/z  - (1/z) P(1/z) / Q(1/z).
		 * There are two rational functions P/Q, one for 0 < y < exp(-32)
		 * and the other for y up to exp(-2).  For larger arguments,
		 * w = y - 0.5, and  x/sqrt(2pi) = w + w**3 R(w**2)/S(w**2)).
		 *
		 *
		 * ACCURACY:
		 *
		 *                      Relative error:
		 * arithmetic   domain        # trials      peak         rms
		 *    DEC      0.125, 1         5500       9.5e-17     2.1e-17
		 *    DEC      6e-39, 0.135     3500       5.7e-17     1.3e-17
		 *    IEEE     0.125, 1        20000       7.2e-16     1.3e-16
		 *    IEEE     3e-308, 0.135   50000       4.6e-16     9.8e-17
		 *
		 *
		 * ERROR MESSAGES:
		 *
		 *   message         condition    value returned
		 * ndtri domain       x <= 0        -MAXNUM
		 * ndtri domain       x >= 1         MAXNUM
		 *
		 */

		/*
		   Cephes Math Library Release 2.8:  June, 2000
		   Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier
		 */


		static double s2pi = 2.50662827463100050242E0;

		static double P0[5] = {
			-5.99633501014107895267E1,
			9.80010754185999661536E1,
			-5.66762857469070293439E1,
			1.39312609387279679503E1,
			-1.23916583867381258016E0,
		};
		static double Q0[8] = {
			1.95448858338141759834E0,
			4.67627912898881538453E0,
			8.63602421390890590575E1,
			-2.25462687854119370527E2,
			2.00260212380060660359E2,
			-8.20372256168333339912E1,
			1.59056225126211695515E1,
			-1.18331621121330003142E0,
		};
		static double P1[9] = {
			4.05544892305962419923E0,
			3.15251094599893866154E1,
			5.71628192246421288162E1,
			4.40805073893200834700E1,
			1.46849561928858024014E1,
			2.18663306850790267539E0,
			-1.40256079171354495875E-1,
			-3.50424626827848203418E-2,
			-8.57456785154685413611E-4,
		};
		static double Q1[8] = {
			1.57799883256466749731E1,
			4.53907635128879210584E1,
			4.13172038254672030440E1,
			1.50425385692907503408E1,
			2.50464946208309415979E0,
			-1.42182922854787788574E-1,
			-3.80806407691578277194E-2,
			-9.33259480895457427372E-4,
		};
		static double P2[9] = {
			3.23774891776946035970E0,
			6.91522889068984211695E0,
			3.93881025292474443415E0,
			1.33303460815807542389E0,
			2.01485389549179081538E-1,
			1.23716634817820021358E-2,
			3.01581553508235416007E-4,
			2.65806974686737550832E-6,
			6.23974539184983293730E-9,
		};
		static double Q2[8] = {
			6.02427039364742014255E0,
			3.67983563856160859403E0,
			1.37702099489081330271E0,
			2.16236993594496635890E-1,
			1.34204006088543189037E-2,
			3.28014464682127739104E-4,
			2.89247864745380683936E-6,
			6.79019408009981274425E-9,
		};
		double ndtri( double y0 )
		{
			double x, y, z, y2, x0, x1;
			int code;
			if( y0 <= 0.0 )
				return( - std::numeric_limits<double>::infinity() );
			if( y0 >= 1.0 )
				return( + std::numeric_limits<double>::infinity() );
			code = 1;
			y = y0;
			if( y > (1.0 - 0.13533528323661269189) ) 
			{
				y = 1.0 - y;
				code = 0;
			}
			if( y > 0.13533528323661269189 )
			{
				y = y - 0.5;
				y2 = y * y;
				x = y + y * (y2 * Polynomialeval( y2, P0, 4)/ Polynomial1eval( y2, Q0, 8 ));
				x = x * s2pi; 
				return(x);
			}
			x = std::sqrt( -2.0 * std::log(y) );
			x0 = x - std::log(x)/x;
			z = 1.0/x;
			if( x < 8.0 ) 
				x1 = z * Polynomialeval( z, P1, 8 )/ Polynomial1eval ( z, Q1, 8 );
			else
				x1 = z * Polynomialeval( z, P2, 8 )/ Polynomial1eval( z, Q2, 8 );
			x = x0 - x1;
			if( code != 0 )
				x = -x;
			return( x );
		}
	} // end namespace Cephes

double normal_quantile(double z, double sigma) {
   // use cephes ndtri function
   return  sigma * Cephes::ndtri(z);
}

//see http://root.cern.ch/root/html/src/RooStats__SamplingDistribution.cxx.html
double InverseCDF(
		vector<double> v,
		double alpha,
		double pvalue, 
		double sigmaVariation, 
		double& inverseWithVariation)
{
	// returns the inverse of the cumulative distribution function, with variations depending on number of samples

	// will need to deal with weights, but for now:
	std::sort(v.begin(), v.end());


	// Acceptance regions are meant to be inclusive of (1-\alpha) of the probability
	// so the returned values of the CDF should make this easy.
	// in particular:
	//   if finding the critical value for a lower bound
	//     when p_i < p < p_j, one should return the value associated with i
	//     if i=0, then one should return -infinity
	//   if finding the critical value for an upper bound
	//     when p_i < p < p_j, one should return the value associated with j
	//     if i = size-1, then one should return +infinity
	//   use pvalue < 0.5 to indicate a lower bound is requested

	// casting will round down, eg. give i
	int nominal = (unsigned int) (pvalue*v.size());

	if(nominal <= 0) {
		inverseWithVariation = -1e20;
		return -1e20;
	}
	else if(nominal >= (int)v.size()-1 ) {
		inverseWithVariation = 1e20;
		return 1e20;
	}
	else if(pvalue < 0.5){
		int delta = (int)(sigmaVariation*sqrt(1.0*nominal)); // note sqrt(small fraction)
		int variation = nominal+delta;

		if(variation>=(int)v.size()-1)
			inverseWithVariation = 1e20;
		else if(variation<=0)
			inverseWithVariation = -1e20;
		else 
			inverseWithVariation =  v[ variation ];

		return v[nominal];
	}
	else if(pvalue >= 0.5){
		int delta = (int)(sigmaVariation*sqrt(1.0*v.size()- nominal)); // note sqrt(small fraction)
		int variation = nominal+delta;

		if(variation>=(int)v.size()-1)
			inverseWithVariation = 1e20;

		else if(variation<=0)
			inverseWithVariation = -1e20;
		else 
			inverseWithVariation =  v[ variation+1 ];

		return v[nominal+1];
	}
	else{
		std::cout << "problem in InverseCDF" << std::endl;
	}
	inverseWithVariation = 1e20;
	return 1e20;

}
};
