#include "UnbinnedBayesianLimit.h"
#include <iostream>
#include <math.h>
using namespace std;
namespace lands{
UnbinnedBayesianLimit::UnbinnedBayesianLimit(){
	debug=0;
	pdfRdmB	= new PdfRandom();
	pdfRdmS	= new PdfRandom();
	_npx=1000;
	_rdm=0;
	_xstop=0;_xstart=0;
	_ns=0; _nb=0;
	_pdfs=0; _pdfb=0; _pars=0; _parb=0;
}
UnbinnedBayesianLimit::~UnbinnedBayesianLimit(){
	delete pdfRdmS;
	delete pdfRdmB;
	_rdm=0;
}
double UnbinnedBayesianLimit::P_vm_Given_brs(double r){
	double retval=exp(-r*_ns);
	double x[1];
	for(int i=0;i<(int)_vmass.size(); i++){
		if(_nb==0) return retval;
		x[0]=_vmass[i];
		double tmp=_pdfb(x,_parb);
		if(tmp==0) continue;
		retval*=(1+r*_ns*_pdfs(x,_pars)/(_nb*tmp));
	}
	return retval;	
}
double UnbinnedBayesianLimit::Integral2(double a, double b, double fEpsilon){
	// - - - from ROOT

	//double fEpsilon = 1.e-8; // configurable  

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
		f1 = P_vm_Given_brs(xx[0]);//(*fcn)(xx, par);


		if (fgAbsValue) {if(f1<0)f1 = -f1;}
		xx[0] = c1-u;
		//f2    = (*fFunction) (xx);
		// --Mingshui
		f2 = P_vm_Given_brs(xx[0]);//(*fcn)(xx, par);

		if (fgAbsValue) {if(f2<0)f2 = -f2;}
		s8   += w[i]*(f1 + f2);
	}
	s16 = 0;
	for (i=4;i<12;i++) {
		u     = c2*x[i];
		xx[0] = c1+u;
		//f1    = (*fFunction) (xx);
		// --Mingshui
		f1 = P_vm_Given_brs(xx[0]);// (*fcn)(xx, par);

		if (fgAbsValue) {if(f1<0)f1 = -f1;}
		xx[0] = c1-u;
		//f2    = (*fFunction) (xx);
		// --Mingshui
		f2 = P_vm_Given_brs(xx[0]);//(*fcn)(xx, par);

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
double UnbinnedBayesianLimit::RLimit(double alpha, double precision, double MinLikelihood, double integralPrecision, double rUpperBound){
	if(_pdfs==0 || _pdfb==0) {
		cout<<"Error: probably you haven't set the the PDF of signal or background"<<endl;
		return 0;
	}
	if(_ns==0 || _nb==0) {
		cout<<"Error: probably you haven't set the total events for signal or background. _ns="<<_ns<<" _nb="<<_nb<<endl;
		return 0;
	}

	double cl_alpha=1-alpha;
	double rmax = rUpperBound;  // r = sigma / sigma_SM
	double norm = Integral2(0, rUpperBound, integralPrecision); // approximately 0~inf

	if(debug) cout<<"UnbinnedBayesianLimit::RLimit norm of likelihood of r ="<<norm<<endl;

	double x[1];
	x[0]=rmax;
	double delta = (P_vm_Given_brs(x[0])/norm-MinLikelihood)/MinLikelihood;


	if(debug) cout<<"UnbinnedBayesianLimit::RLimit delta="<<delta<<endl;
	while(delta>0){

		rmax=2*rmax;
		x[0]=rmax;
		delta = (P_vm_Given_brs(x[0])/norm-MinLikelihood)/MinLikelihood;
		if(debug) cout<<"UnbinnedBayesianLimit::RLimit trying to let delta <=0, current delta="<<delta<<"  x[0]="<<x[0]<<endl;
	}

	double x1=0, x2=rmax;
	while (fabs(delta) > precision) {
		if(delta>0) x1=x[0];
		else x2=x[0];
		x[0]=(x1+x2)/2.;
		delta = (P_vm_Given_brs(x[0])/norm-MinLikelihood)/MinLikelihood;	  
		if(debug) cout<<"UnbinnedBayesianLimit::RLimit trying to let |delta|<="<<precision<<", current |delta|="<<fabs(delta)<<"  x[0]="<<x[0]<<endl;
	}
	rmax = x[0];

	if(debug) cout<<"UnbinnedBayesianLimit::RLimit rmax="<<rmax<<endl;

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


	//	double norm = Integral2(0.,rmax);
	double r1 = 0;
	double r2 = rmax;
	rmax = (r1 + r2)/2.;
	delta = Integral2(0,rmax, integralPrecision)/norm - cl_alpha;
	while (fabs(delta) > precision) 
	{
		if (delta < 0) r1 = rmax;
		else r2 = rmax;
		rmax = (r1 + r2)/2.;
		delta = Integral2(0,rmax, integralPrecision)/norm - cl_alpha;
		if(debug) cout<<"UnbinnedBayesianLimit::RLimit trying to let |integral - alpha|<="<<precision
			<<", current |delta|="<<fabs(delta)<<"  integral from "<<0<<" to "<<rmax<<endl;
	} 
	if(debug) cout<<"UnbinnedBayesianLimit::RLimit final rmax="<<rmax<<endl;
	_r95=rmax;
	return _r95;
}
void UnbinnedBayesianLimit::SetPdfs(double (*fcn)(double *, double *), double *par){
	_pdfs=fcn;
	_pars=par;
	if(debug){
		double x[1];
		x[0]=130;
		cout<<"UnbinnedBayesianLimit::SetPdfs  testing pdfs("<<x[0]<<")="<<_pdfs(x,_pars)<<endl;
	}
}
void UnbinnedBayesianLimit::SetPdfb(double (*fcn)(double *, double *), double *par){
	_pdfb=fcn;
	_parb=par;
	if(debug){
		double x[1];
		x[0]=130;
		cout<<"UnbinnedBayesianLimit::SetPdfb  testing pdfb("<<x[0]<<")="<<_pdfb(x,_parb)<<endl;
	}
}
void UnbinnedBayesianLimit::SetVmass(double *v, int nevts){
	_vmass.clear();
	for(int i=0; i<nevts; i++){
		_vmass.push_back(v[i]);
	}
}
void UnbinnedBayesianLimit::SetVmass(vector<double> v){
	_vmass.clear();
	for(int i=0; i<(int)v.size(); i++){
		_vmass.push_back(v[i]);
	}
}
void UnbinnedBayesianLimit::SetExpectedNsNb(double ns, double nb){
	_ns=ns;
	_nb=nb;
}
double UnbinnedBayesianLimit::R_3A(double alpha, double precision, double MinLikelihood, double rUpperBoud){
	if(_ns==0 || _nb==0) {
		cout<<"Error: probably you haven't set the total events for signal or background. _ns="<<_ns<<" _nb="<<_nb<<endl;
		return 0;
	}
	if(_xstop==_xstart) {
		cout<<"Error: probably you haven't set the start and stop values to be tested. _xstart="<<_xstart<<" _xstop"<<_xstop<<endl;
		return 0;
	}
	if(_pdfs==0 || _pdfb==0) {
		cout<<"Error: probably you haven't set the the PDF of signal or background"<<endl;
		return 0;
	}

	// --- part 3A
	double par[3000];
	int nbins=1000;  //need to be configurable
	double dm = (_xstop-_xstart)/(double)nbins;
	for(int i=0; i<nbins; i++){
		// --- calculate delta area of dm by integral 
		par[0+i*3]=Integral(_pdfs, _xstart+i*dm, _xstart+(i+1)*dm, _pars)*_ns;
		par[1+i*3]=Integral(_pdfb, _xstart+i*dm, _xstart+(i+1)*dm, _parb)*_nb;
		par[2+i*3]=par[1+i*3];
	}
	double limit3A = R_CL(par, 1-alpha, nbins, precision, MinLikelihood, rUpperBoud);      // 1 channel
	return limit3A;
}
void UnbinnedBayesianLimit::RunMCexps(int nexps, double alpha, double precision, double MinLikelihood, double integralPrecision, double rUpperBoud){
	if(_ns==0 || _nb==0) {
		cout<<"Error: probably you haven't set the total events for signal or background. _ns="<<_ns<<" _nb="<<_nb<<endl;
		return ;
	}
	if(_xstop==_xstart) {
		cout<<"Error: probably you haven't set the start and stop values to be tested. _xstart="<<_xstart<<" _xstop"<<_xstop<<endl;
		return ;
	}
	if(_pdfs==0 || _pdfb==0) {
		cout<<"Error: probably you haven't set the the PDF of signal or background"<<endl;
		return ;
	}
	if(!_rdm) {
		cout<<"Error: probably you haven't set the random engine"<<endl;
		return ;
	}
	pdfRdmS->SetFunction(_pdfs, _pars);
	pdfRdmS->SetRndGen(_rdm);
	pdfRdmS->SetNpx(_npx);
	pdfRdmS->SetRange(_xstart, _xstop);
	pdfRdmB->SetFunction(_pdfb, _parb);
	pdfRdmB->SetRndGen(_rdm);
	pdfRdmB->SetNpx(_npx);
	pdfRdmB->SetRange(_xstart, _xstop);

	double r;
	_vR.clear();
	if(debug)cout<<"UnbinnedBayesianLimit::RunMCexps total "<<nexps<<" exps"<<endl;
	for(int i=0; i<nexps; i++){
		if(debug)cout<<" UnbinnedBayesianLimit::RunMCexps "<<i<<"th exp ......"<<endl;
		int ni = _rdm->Poisson(_nb);
		_vmass.clear();
		for(int j=0; j<ni; j++){
			_vmass.push_back(pdfRdmB->GetRandom());		
		}
		r=RLimit(alpha, precision, MinLikelihood, integralPrecision, rUpperBoud); //precision of integration is not enough
		_vR.push_back(r);
		if(debug) {
			cout<<"limit in "<<i<<"th exp is "<<r<<", throwing "<<ni<<" events, they are: ";
			for(int j=0; j<ni; j++) cout<<_vmass[j]<<" ";
			cout<<endl;
		}	
	}	
}
double UnbinnedBayesianLimit::Get_R_mean(){
	double r=0;
	int nexps=_vR.size();
	if(nexps==0) {cout<<"UnbinnedBayesianLimit::Get_R_mean  Warning: no pseudoexps yet"<<endl; return 0;}
	for(int i=0; i<nexps; i++ ){
		r+=_vR[i];
	}
	r=r/(double)nexps;
	return r;	
}
double UnbinnedBayesianLimit::Get_R_meanerr(){
	double r=Get_R_mean();
	double rs=0, rerr=0; 
	int nexps=_vR.size();
	if(nexps==0) {cout<<"UnbinnedBayesianLimit::Get_R_mean  Warning: no pseudoexps yet"<<endl; return 0;}
	for(int i=0; i<nexps; i++ ){
		rs+=(_vR[i]-r)*(_vR[i]-r);
	}
	rerr=sqrt(rs)/(double)nexps;
	return rerr;	
}
};
