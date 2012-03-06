#include "MLLxsBands.h"
#include <time.h> // upto micro second
#include <iostream>
#include <cstdio>
#include "TMath.h"
using std::cout;
using std::endl;
namespace lands{
	MLLxsBands::MLLxsBands(CountingModel* cms){
		//_frequentist=cb;
		_cms=cms;
		_debug=0; 
	}
	MLLxsBands::~MLLxsBands(){}
	void MLLxsBands::Bands(int noutcome){
		_noutcomes=noutcome;
		Bands();
	}
	void MLLxsBands::Bands(){
		clock_t start_time=clock(), cur_time=clock(), funcStart_time=clock();
		if(!_cms) {cout<<" Bands, model not set yet"<<endl; return; }
		//if(!_frequentist)  {cout<<" Bands, frequentist not set yet"<<endl; return; }		


		double rmeancls=0;

		_vrcls.clear(); _vpcls.clear();// _difrcls.clear();

		vector<double> vrcls; vrcls.clear();
		vector<double> vpcls; vpcls.clear();

		cout<<" MLLxs = "<<_cms->Get_fittedParsInData_sb()[0]<<endl;
		_cms->SetSignalScaleFactor(_cms->Get_fittedParsInData_sb()[0]);

		double *pars = new double[_cms->Get_max_uncorrelation()+1];
		for(int i=0; i<_noutcomes; i++){
			vdata_global = (VDChannel)_cms->GetToyData_H1(_cms->Get_fittedParsInData_sb());
			if(_cms->hasParametricShape()){
				_cms->SetTmpDataForUnbinned(_cms->Get_v_pdfs_roodataset_toy());
			}
			double tmpr=0, tmperr=0; 
			double ErrorDef = TMath::ChisquareQuantile(0.68 , 1);// (confidenceLevel, ndf)
			int success[1]={0};

			double y0_2 =  MinuitFit(102, tmpr, tmperr, ErrorDef, pars, false, _debug, success) ;  //202, 201, 2:  allow mu to be negative.   102  constrain to be positive 
			if(_debug)cout<<"Toy#"<<i<<": L="<<y0_2<<" muhat="<<pars[0]<<", from minos fit asymmetry 68% CL:  ["<<tmperr<<","<<tmpr<<"]"<<endl; // upperL, lowerL
			double mu_hat = pars[0];

			vrcls.push_back(mu_hat);
			vpcls.push_back(1./_noutcomes);
			rmeancls+=mu_hat/_noutcomes;

		}


		if(1) {
			int nr = vrcls.size(); 
			double *r_v = new double[nr];	
			int *ir_v = new int[nr];	
			for(int i=0; i<nr; i++) r_v[i]=vrcls.at(i);
			Sort(nr, r_v, ir_v, 0); // from small to large
			for(int i=0; i<nr; i++){
				_vrcls.push_back(vrcls.at(ir_v[i]));
				if(i==0) _vpcls.push_back(vpcls.at(ir_v[i]));
				else _vpcls.push_back( _vpcls.at(i-1) + vpcls.at(ir_v[i]) ); 
			}	
			if(r_v) delete [] r_v;
			if(ir_v) delete [] ir_v;

			if(_debug) {
				cout<<" ProjectingMLLxs, printing out the muhat and cummulative p for all possible outcomes:"<<endl;
				for(int i=0; i<(int)_vrcls.size(); i++)
					printf("r %5.2f p %5.4f\n",_vrcls[i], _vpcls[i]);
			}

			GetBandsByFermiCurveInterpolation(_vrcls,_vpcls, _rcls[1], _rcls[3], _rcls[0], _rcls[4]);
			_rcls[5]=rmeancls; _rcls[2]=GetBandByFermiCurveInterpolation(_vrcls, _vpcls, 0.5);
			if(_debug) {
				cout<<"\t ProjectingMLLxss projecting muhat at "<<endl; 
				cout<<"    -2 sigma ="<<_rcls[0]<<endl;
				cout<<"    -1 sigma ="<<_rcls[1]<<endl;
				cout<<"    median   ="<<_rcls[2]<<endl;
				cout<<"     1 sigma ="<<_rcls[3]<<endl;
				cout<<"     2 sigma ="<<_rcls[4]<<endl;
				cout<<"     mean    ="<<_rcls[5]<<endl;
			}
		}
		if(_debug) { start_time=cur_time; cur_time=clock(); cout << "\t\t\t TIME_in_BAND, final printing took " <<(cur_time - start_time)/1000. << " milisec\n"; }
		if(_debug) { cur_time=clock(); cout << "\t\t\t TIME_in_BAND finishing "<<_noutcomes<<" outcomes: " << (cur_time - funcStart_time)/1000000./60. << " minutes\n"; }
	}	
};
