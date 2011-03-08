#include "SignificanceBands.h"
#include <time.h> // upto micro second
#include <iostream>
#include <cstdio>
using std::cout;
using std::endl;
namespace lands{
	SignificanceBands::SignificanceBands(CLsBase *cb, CountingModel* cms){
		_frequentist=cb; _cms=cms;
		_debug=0; 
	}
	SignificanceBands::~SignificanceBands(){}
	void SignificanceBands::Bands(int noutcome){
		_noutcomes=noutcome;
		Bands();
	}
	void SignificanceBands::Bands(){
		clock_t start_time=clock(), cur_time=clock(), funcStart_time=clock();
		if(!_cms) {cout<<" Bands, model not set yet"<<endl; return; }
		if(!_frequentist)  {cout<<" Bands, frequentist not set yet"<<endl; return; }		

		_cms->SetSignalScaleFactor(1.);

		double rmeancls=0;

		_vrcls.clear(); _vpcls.clear(); _difrcls.clear();

		vector<double> vrcls; vrcls.clear();
		vector<double> vpcls; vpcls.clear();

		if(_debug >= 10){ // print out s,b se, be, correlation, err_pdf_type
			cout<<"\t DoingStatisticalBandsForSignificance print out details on signal, background and uncertainties ...."<<endl;
			_cms->Print();
		}

		vector< VIChannel > vvPossibleOutcomes; vvPossibleOutcomes.clear();
		vector<int> vEntries; vEntries.clear();  // index of vvPossibleOutcomes are the same as vEntries
		VIChannel vbkg_tmp; 

		if(_debug) cout<<" _noutcomes="<<_noutcomes<<endl;

		for(int n=0; n<_noutcomes; n++){
			vbkg_tmp.clear(); 
			vbkg_tmp=_cms->GetToyData_H0(); //// this is important, need throw H1 hypothesis
			if(_debug>=10) {
				cout<<"H0 : ";
				for(int j=0; j<vbkg_tmp.size(); j++){
					if(_debug>=10)	printf("%4d ", vbkg_tmp.at(j));
				}
				cout<<endl;		
			}
			vbkg_tmp=_cms->GetToyData_H1(); //// this is important, need throw H1 hypothesis
			if(_debug>=10) {
				cout<<"H1 : ";
				for(int j=0; j<vbkg_tmp.size(); j++){
					if(_debug>=10)	printf("%4d ", vbkg_tmp.at(j));
				}
				cout<<endl;		
			}
			bool flag=true;
			for(int t=0; t<vvPossibleOutcomes.size(); t++){
				if(vbkg_tmp==vvPossibleOutcomes.at(t)) {
					flag=false;
					vEntries.at(t)+=1;
					break;
				}
			}
			if(flag) { 
				vvPossibleOutcomes.push_back(vbkg_tmp);
				vEntries.push_back(1);	
			}
		}	// _noutcomes

		if(_debug) { start_time=cur_time; cur_time=clock(); cout << "\t\t\t TIME_in_BAND, generating "<<_noutcomes<<" outcomes: " << (cur_time - start_time)/1000000. << " sec\n"; }

		// rank the possible outcomes and sort it
		int npsbl_outcomes = vvPossibleOutcomes.size();
		int *Entries_v = new int[npsbl_outcomes];
		for(int i=0; i<npsbl_outcomes; i++){
			Entries_v[i]=vEntries.at(i);
		}
		int *iEntries_v = new int[npsbl_outcomes];
		Sort(npsbl_outcomes, Entries_v, iEntries_v, 1); // sorted it by "down",  Desc
		if(_debug) {
			cout<<"\n\t\t ProjectingLimits possible outcomes and their entries "<<endl;
			cout<<"\t\t n_possible_outcomes size= "<<npsbl_outcomes<<endl;
			int step = npsbl_outcomes/10;
			for(int i=0; i<npsbl_outcomes; i+=(step+1) ){
				printf("\t%dth\t ", i);
				for(int j=0; j<vvPossibleOutcomes.at(iEntries_v[i]).size(); j++){
					if(_debug>=10)	printf("%4d ", vvPossibleOutcomes.at(iEntries_v[i]).at(j));
				}		
				printf("\t %5d\n", vEntries.at(iEntries_v[i]));
			}
		}

		double nchs=_cms->NumOfChannels();
		int tmpcount;
		if(npsbl_outcomes<50) tmpcount = int(npsbl_outcomes/10);
		else {
			tmpcount = int(npsbl_outcomes/100);
		}

		cout<<"\t Computing expected Upper Limit for "<<_noutcomes<<" sets of outcomes "<<endl;
		double r, tmp;
		for(int n=0; n<npsbl_outcomes; n++){
			if( tmpcount==0 || (tmpcount!=0 && (n%tmpcount == 0)) ) printf(" ... ... ... process %3.0f\%\n", n/(double)npsbl_outcomes*100);
			double p=vEntries.at(iEntries_v[n])/(double)_noutcomes;
			if(_debug){
				cout<<"\t ProjectingRLimits scanning outcomes: "<<n<<"th \t";
				for(int j=0; j<nchs; j++){
					if(_debug>=10)printf("%4d ", vvPossibleOutcomes.at(iEntries_v[n]).at(j));
				}		
				cout<<"\t p="<<p<<endl;
			}
			//_cms->SetData( vvPossibleOutcomes.at(iEntries_v[n]) ); 
			cms_global = _cms;
			vdata_global = vvPossibleOutcomes.at(iEntries_v[n]);

			//_cms->SetData(vdata_global);
			//_cms->Print();

			int nentries_for_thisR=(int)(vEntries.at(iEntries_v[n]));
			//int nentries_for_thisR=(int)(vEntries.at(iEntries_v[n]));
			if(1) { 
				//--------------calc r95% with (s,b,n) by throwing pseudo experiments and using the fits to get r95% at CLs=5
				double m2lnQ = MinuitFit(3,tmp, tmp) - MinuitFit(2, tmp, tmp);
				if(m2lnQ < 0) cout<<"m2lnQ_profiled < 0, ="<<m2lnQ<<endl;
				r = sqrt(fabs(m2lnQ));
				vrcls.push_back(r);
				vpcls.push_back(p);
				rmeancls+=r*p;
				for(int ntmp=0; ntmp<nentries_for_thisR; ntmp++){
					_difrcls.push_back(r);
				}
				if(_debug)cout<<"n="<<n<<" cls_r= "<<r<<" p= "<<p<<" repeated "<<nentries_for_thisR<<endl;
				if(_debug) { start_time=cur_time; cur_time=clock(); cout << "\t\t\t TIME_in_BAND, doCLs "<< n << " took " <<(cur_time - start_time)/1000000. << " sec\n"; }
			}
			fflush(stdout);
		}

		if(Entries_v) delete [] Entries_v;
		if(iEntries_v) delete [] iEntries_v;

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
				cout<<" ProjectingSignificance, printing out for cls the r and cummulative p for all possible outcomes:"<<endl;
				for(int i=0; i<(int)_vrcls.size(); i++)
					printf("r %5.2f p %5.4f\n",_vrcls[i], _vpcls[i]);
			}

			GetBandsByFermiCurveInterpolation(_vrcls,_vpcls, _rcls[1], _rcls[3], _rcls[0], _rcls[4]);
			_rcls[5]=rmeancls; _rcls[2]=GetBandByFermiCurveInterpolation(_vrcls, _vpcls, 0.5);
			if(_debug) {
				cout<<"\t ProjectingSignificances projecting cls r at "<<endl; 
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
