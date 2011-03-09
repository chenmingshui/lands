#include "UnbinnedCLsLimit.h"
#include "UnbinnedBayesianLimit.h"
#include <iostream>
#include <math.h>
using namespace std;
namespace lands{
UnbinnedCLsLimit::UnbinnedCLsLimit(){
	debug=0;
	pdfRdmB	= new PdfRandom();
	pdfRdmS	= new PdfRandom();
	_cls = new CLsBase();
	_npx=1000;

	_rdm=0;
	_xstop=0;_xstart=0;
	_ns=0; _nb=0;
	_pdfs=0; _pdfb=0; _pars=0; _parb=0;
}
UnbinnedCLsLimit::~UnbinnedCLsLimit(){
	delete pdfRdmS;
	delete pdfRdmB;
	delete _cls;
	_rdm=0;
}
double UnbinnedCLsLimit::lnQ(vector<double> vmass){
	double retval=0;
	double x[1];
	for(int i=0;i<(int)vmass.size(); i++){
		
		if(_nb==0) return retval;
		x[0]=vmass[i];
		double tmp=_pdfb(x,_parb);
		if(tmp==0){
			cout<<"Error: pdfb("<<x<<")=0"<<endl;	
			continue;
		}
		tmp=(1+_ns*_pdfs(x,_pars)/(_nb*tmp));
		if(tmp<=1) {
			cout<<"Error: 1+r*Stot*pdfs(x)/(Btot*pdfb(x))="<<tmp<<"<=1"<<endl;
			continue;
		}
		retval+=log(tmp);
	}
	retval-=_ns;

//	if(debug) cout<<"UnbinnedCLsLimit::lnQ lnQ="<<retval<<endl;
	return retval;	
}
double UnbinnedCLsLimit::CLs(int nexps){
	// before running, should check all is ok....... such as _rdm is ok, _xstart and _xstop is set  ..... ..... ...
	if(CheckOk()==false) return 0;
	double retval=0;

	if(debug) cout<<"UnbinnedCLsLimit::CLs starting.."<<endl;	
		
	double lnQ_data=lnQ(_vmass);

	if(debug) cout<<"UnbinnedCLsLimit::CLs lnQ_data="<<lnQ_data<<endl;	
	vector<double> vlnQ_b, vlnQ_sb;
	vlnQ_b.clear();vlnQ_sb.clear();
	int quarter = nexps/4;
	for(int i=0; i<nexps; i++){
		vector<double> vmassi;
		vmassi.clear();
		int ni=_rdm->Poisson(_nb);
		for(int j=0; j<ni; j++){
			vmassi.push_back(pdfRdmB->GetRandom());		
		}
		vlnQ_b.push_back(lnQ(vmassi));		
	
		vmassi.clear();
		ni=_rdm->Poisson(_nb);
		for(int j=0; j<ni; j++){
			vmassi.push_back(pdfRdmB->GetRandom());		
		}
		ni=_rdm->Poisson(_ns);
		for(int j=0; j<ni; j++){
			vmassi.push_back(pdfRdmS->GetRandom());		
		}
		vlnQ_sb.push_back(lnQ(vmassi));		

		if(debug && i%quarter==1) cout<<" UnbinnedCLsLimit::CLs running "<<i<<"th experiment"<<endl;	
	}
	
	_cls->SetDebug(debug);
	_cls->SetLogQ_b(vlnQ_b);
	_cls->SetLogQ_sb(vlnQ_sb);
	_cls->SetLogQ_data(lnQ_data);
	double err;
	retval=_cls->CLs(err);		

	if(debug) cout<<"UnbinnedCLsLimit::CLs done, CLs="<<retval<<" +/- "<<err<<endl;	

	return retval;
}
double UnbinnedCLsLimit::RLimit(double alpha, double epsilon, int nexps){
	if(CheckOk()==false) return 0;
	double _r95=0;
	double ns_dontChange=_ns;
	_vR.clear(); _vCLs.clear();

	cout<<endl<<"UnbinnedCLsLimit::RLimit looking for C.L. 95% Limit on the ratio ----"<<endl;
	UnbinnedBayesianLimit ubl;
	ubl.SetDebug(debug);
	ubl.SetPdfs(_pdfs, _pars);
	ubl.SetPdfb(_pdfb, _parb);
	ubl.SetExpectedNsNb(_ns, _nb);
	ubl.SetVmass(_vmass);
	double r0,r1;
	r0 = ubl.RLimit(alpha);  
	cout<<"UnbinnedCLsLimit::RLimit start with estimate from Bayesian technique, r95%="<<r0<<endl;

	_ns=r0*ns_dontChange;
	double cl0=CLs(nexps);
	_vR.push_back(r0);_vCLs.push_back(cl0);
	cout<<"r="<<r0<<"  alpha="<<cl0<<endl;
	if(fabs(cl0-alpha)<=epsilon) {
		_r95=r0;
		cout<<"Converge at alpha="<<alpha<<"+/-"<<epsilon<<" by "<<"1 iteration"<<endl;
		_ns=ns_dontChange;
		return _r95;
	}

	r1=r0*0.90; //----------usually, CLs-limit is more aggresive than Bayesian's, about 10% smaller.
	_ns=r1*ns_dontChange;
	double cl1=CLs(nexps);
	_vR.push_back(r1);_vCLs.push_back(cl1);
	cout<<"r="<<r1<<"  alpha="<<cl1<<endl;
	if(fabs(cl1-alpha)<=epsilon) {
		_r95=r1;
		cout<<"Converge at alpha="<<alpha<<"+/-"<<epsilon<<" by "<<"2 iterations"<<endl;
		_ns=ns_dontChange;
		return _r95;
	}

	bool foundit=false;
	double rmid=0;
	int nmaxrepeat=0;	
	while(!foundit && nmaxrepeat<=30 ){ 
		// --- -  --
		//  using linear interpolation to do converge will be quicker
		// ---------
		rmid=LinearInterpolation(r0,cl0,r1,cl1,alpha);
		//	rmid= (alpha-(cl1-((cl1-cl0)/(r1-r0))*r1))/((cl1-cl0)/(r1-r0));
		if(debug)cout<<" r0="<<r0<<" cl0="<<cl0<<" r1="<<r1<<" cl1="<<cl1<<" rmid="<<rmid<<endl;
		_ns=rmid*ns_dontChange;
		double clmid=CLs(nexps);
		_vR.push_back(rmid);_vCLs.push_back(clmid);
		cout<<"r="<<rmid<<"  alpha="<<clmid<<endl;
		if(fabs(clmid-alpha)<epsilon) foundit=true; //alpha=0.05 C.L. 95% 		
		else {
			double x[3]={r0,rmid,r1}; double y[3]={cl0,clmid,cl1};		
			int iy[3];
			Sort(3,y,iy,0);
			if(alpha<y[iy[1]]){ //---------kick out a number among r0, r1, rmid
				cl0=y[iy[0]]; cl1=y[iy[1]];
				r0=x[iy[0]]; r1=x[iy[1]];
			}
			else if(alpha>y[iy[1]]){
				cl0=y[iy[1]]; cl1=y[iy[2]];
				r0=x[iy[1]]; r1=x[iy[2]];
			}
			else {
				cout<<"CANNOT converge: r0="<<r0<<" cl0="<<cl0<<" r1="<<r1<<" cl1="<<cl1<<" rmid="<<rmid<<" clmid="<<clmid<<endl;
			}
		}
		nmaxrepeat++;
		//------------------------????????????? is CLs mono-increased/decreased as r ?  has to investigate it 
		if(!foundit){
			if(rmid==_vR[_vR.size()-1]) foundit=true; 
			cout<<"We get rmid="<<rmid<<" twice, and decide to stop here...."<<endl;
		}
		if(foundit && _vR.size()>1){	
			bool hasCLsGT05=false;
			bool hasCLsLT05=false;
			for(int icls= 0; icls<_vR.size(); icls++){
				if(_vCLs[icls]<=alpha) hasCLsLT05=true;
				if(_vCLs[icls]>=alpha) hasCLsGT05=true;
			}	
			if(!hasCLsLT05) {
				rmid= rmid*1.05;
				_ns=rmid*ns_dontChange;
				double clmid=CLs(nexps);
				_vR.push_back(rmid);_vCLs.push_back(clmid);
				cout<<"r="<<rmid<<"  alpha="<<clmid<<endl;
			}
			if(!hasCLsGT05) {
				rmid= rmid*0.95;
				_ns=rmid*ns_dontChange;
				double clmid=CLs(nexps);
				_vR.push_back(rmid);_vCLs.push_back(clmid);
				cout<<"r="<<rmid<<"  alpha="<<clmid<<endl;
			}
		}
	}
	int nsize=_vR.size();

	if(!foundit){cout<<"Not converge to "<<alpha<<"+/-"<<epsilon<<" with "<<nsize<<" iterations"<<endl;}
	else cout<<"Converge at alpha="<<alpha<<"+/-"<<epsilon<<" by "<<nsize<<" iterations"<<endl;	

	if(nsize==1) {_r95=rmid; _ns=ns_dontChange; return _r95;}

	// --- -get the _r95 @ alpha=0.05 by linear interpolation
	double x1=0, y1=0, x2=0, y2=0;

	double *dCLs=new double[nsize];
	for(int i=0; i<nsize; i++){
		dCLs[i]=_vCLs[i];
	}

	int *ix=new int[nsize];
	Sort(nsize, dCLs, ix, 0);
	for(int i=0; i<nsize; i++){
		if(dCLs[ix[i]]==alpha) {_r95=_vR[ix[i]]; _ns=ns_dontChange; return _r95;}
		if(dCLs[ix[i]]<alpha) { 
			x1=_vR[ix[i]]; 
			y1=dCLs[ix[i]];
		}
		if(dCLs[ix[i]]>alpha && x2==0){
			x2=_vR[ix[i]]; 
			y2=dCLs[ix[i]];
		}
	}
	if(!x1 || !x2 || !y1 || !y2)
		cout<<"ERROR x1="<<x1<<" y1="<<y1<<" x2="<<x2<<" y2="<<y2<<endl;
	_r95=LinearInterpolation(x1,y1,x2,y2,alpha);
	cout<<"final upperlimit on r is "<<_r95<<endl;

	_ns=ns_dontChange;
	return _r95;	
}
double UnbinnedCLsLimit::RLimit(double alpha, double epsilon, double rmin, double rmax, int nexps){
	if(CheckOk()==false) return 0;
	double _r95=0;
	double ns_dontChange=_ns;

	_vR.clear(); _vCLs.clear();

	cout<<endl<<"UnbinnedCLsLimit::RLimit looking for C.L. 95% Limit on the ratio ----"<<endl;

	double r0=rmin, r1=rmax;

	_ns=r0*ns_dontChange;
	double cl0=CLs(nexps);
	_vR.push_back(r0);_vCLs.push_back(cl0);
	cout<<"r="<<r0<<"  alpha="<<cl0<<endl;

	_ns=r1*ns_dontChange;
	double cl1=CLs(nexps);
	_vR.push_back(r1);_vCLs.push_back(cl1);
	cout<<"r="<<r1<<"  alpha="<<cl1<<endl;

	bool foundit=false;
	int nmaxrepeat=0;
	while ( !foundit && nmaxrepeat<30){
		nmaxrepeat+=1;

		// --- -  --
		//  using linear interpolation to do converge will be quicker
		// ---------
		double rmid = r0+0.5*(r1-r0);
		_ns=rmid*ns_dontChange;
		double clmid=CLs(nexps);
		_vR.push_back(rmid);_vCLs.push_back(clmid);
		cout<<"r="<<rmid<<"  alpha="<<clmid<<endl;
		if(fabs(clmid-alpha)<epsilon) foundit=true; //alpha=0.05 C.L. 95% 		
		else if(clmid>alpha && cl1<alpha) {r0=rmid; cl0=clmid;}
		else if(clmid<alpha && cl1<alpha) {r1=rmid; cl1=clmid;}
		else if(clmid<alpha && cl1>alpha) {r0=rmid; cl0=clmid;} // cl1>alpha then wrong ?
		else if(clmid<alpha && cl1>alpha) {r1=rmid; cl1=clmid;} 
	}
	int nsize=_vR.size();
	cout<<"Converge at alpha="<<alpha<<"+/-"<<epsilon<<" by "<<nsize<<" iterations"<<endl;	

	// --- -get the _r95 @ alpha=0.05 by linear interpolation
	double x1=0, y1=0, x2=0, y2=0;

	double *dCLs=new double[nsize];
	for(int i=0; i<nsize; i++){
		dCLs[i]=_vCLs[i];
	}

	int *ix=new int[nsize];
	Sort(nsize, dCLs, ix, 0);
	for(int i=0; i<nsize; i++){
		if(dCLs[ix[i]]==alpha) {_r95=_vR[ix[i]]; _ns=ns_dontChange; return _r95;}
		if(dCLs[ix[i]]<alpha) { 
			x1=_vR[ix[i]]; 
			y1=dCLs[ix[i]];
		}
		if(dCLs[ix[i]]>alpha && x2==0){
			x2=_vR[ix[i]]; 
			y2=dCLs[ix[i]];
		}
	}
	_r95=LinearInterpolation(x1,y1,x2,y2,alpha);

	_ns=ns_dontChange;
	return _r95;	
}
void UnbinnedCLsLimit::SetPdfs(double (*fcn)(double *, double *), double *par){
	_pdfs=fcn;
	_pars=par;
	if(debug){
		double x[1];
		x[0]=130;
		cout<<"UnbinnedCLsLimit::SetPdfs  testing pdfs("<<x[0]<<")="<<_pdfs(x,_pars)<<endl;
	}
}
void UnbinnedCLsLimit::SetPdfb(double (*fcn)(double *, double *), double *par){
	_pdfb=fcn;
	_parb=par;
	if(debug){
		double x[1];
		x[0]=130;
		cout<<"UnbinnedCLsLimit::SetPdfb  testing pdfb("<<x[0]<<")="<<_pdfb(x,_parb)<<endl;
	}
}
void UnbinnedCLsLimit::SetVmass(double *v, int nevts){
	_vmass.clear();
	for(int i=0; i<nevts; i++){
		_vmass.push_back(v[i]);
	}
}
void UnbinnedCLsLimit::SetVmass(vector<double> v){
	_vmass.clear();
	for(int i=0; i<(int)v.size(); i++){
		_vmass.push_back(v[i]);
	}
}
void UnbinnedCLsLimit::SetExpectedNsNb(double ns, double nb){
	_ns=ns;
	_nb=nb;
}
bool UnbinnedCLsLimit::CheckOk(){
	bool ok=true;
	if(_ns==0 || _nb==0) {
		cout<<"Error: probably you haven't set the total events for signal or background. _ns="<<_ns<<" _nb="<<_nb<<endl;
		ok=false;
	}
	if(_xstop==_xstart) {
		cout<<"Error: probably you haven't set the start and stop values to be tested. _xstart="<<_xstart<<" _xstop"<<_xstop<<endl;
		ok=false;
	}
	if(_pdfs==0 || _pdfb==0) {
		cout<<"Error: probably you haven't set the the PDF of signal or background"<<endl;
		ok=false;
	}
	if(!_rdm) {
		cout<<"Error: probably you haven't set the random engine"<<endl;
		ok=false;
	}
	if(ok){
		pdfRdmS->SetFunction(_pdfs, _pars);
		pdfRdmS->SetRndGen(_rdm);
		pdfRdmS->SetNpx(_npx);
		pdfRdmS->SetRange(_xstart, _xstop);
		pdfRdmB->SetFunction(_pdfb, _parb);
		pdfRdmB->SetRndGen(_rdm);
		pdfRdmB->SetNpx(_npx);
		pdfRdmB->SetRange(_xstart, _xstop);
	}
	return ok;
}
};
