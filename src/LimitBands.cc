#include "LimitBands.h"
#include <time.h> // upto micro second
#include <cstdio>
#include <iostream>
#include "RooAbsData.h"
#include "RooDataSet.h"
#include "TTree.h"
#include <map>
#include "UtilsROOT.h"
#include "PlotUtilities.h"

using namespace RooFit;
namespace lands{
	LimitBands::LimitBands(){
		_frequentist=0;	_clslimit=0; _byslimit=0; _cms=0;
		_doCLs=0; _doBys=0; _debug=0; _plotLevel=0;bOnlyEvalCL_forVR = 0;
		bQuickEstimateInitialLimit = false;
		initialRmin = 0.01; initialRmax=1; nSteps=3;
	}
	LimitBands::LimitBands(CLsLimit *cl, CLsBase *cb, CountingModel* cms){
		_frequentist=cb;	_clslimit=cl;  _byslimit=0; _cms=cms;
		_doCLs=0; _doBys=0; _debug=0;  _plotLevel=0; bOnlyEvalCL_forVR = 0;
		bQuickEstimateInitialLimit = false;
		initialRmin = 0.01; initialRmax=1; nSteps=3;
	}
	LimitBands::LimitBands(BayesianBase *bb, CountingModel* cms){
		_frequentist=0;	_clslimit=0; _byslimit=bb; _cms=cms;
		_doCLs=0; _doBys=0; _debug=0;  _plotLevel=0; bOnlyEvalCL_forVR = 0;
		bQuickEstimateInitialLimit = false;
		initialRmin = 0.01; initialRmax=1; nSteps=3;
	}
	LimitBands::LimitBands(CLsLimit *cl, CLsBase *cb, BayesianBase *bb, CountingModel* cms){
		_frequentist=cb;	_clslimit=cl;  _byslimit=bb; _cms=cms; bb->SetModel(cms);
		_doCLs=0; _doBys=0; _debug=0;  _plotLevel=0; bOnlyEvalCL_forVR = 0;
		bQuickEstimateInitialLimit = false;
		initialRmin = 0.01; initialRmax=1; nSteps=3;
	}
	LimitBands::~LimitBands(){}
	void LimitBands::CLsLimitBands(double alpha, int noutcome, int ntoysM2lnQ){
		fAlpha=alpha; fConfidenceLevel=1-fAlpha; _noutcomes=noutcome;
		_doCLs=true; _ntoysM2lnQ=ntoysM2lnQ; 
		if(_clslimit){ _clslimit->SetAlpha(fAlpha); }
		Bands();
	}
	void LimitBands::BysLimitBands(double alpha, int noutcome, int ntoysBys){
		fAlpha=alpha; fConfidenceLevel=1-fAlpha; _noutcomes=noutcome; _actualOutComesForBys = _noutcomes;
		_doBys=true; _ntoysBys=ntoysBys; 
		if(_byslimit){_byslimit->SetAlpha(fAlpha); _byslimit->SetNumToys(_ntoysBys); }
		Bands();
	}
	void LimitBands::Bands(double alpha, int noutcome, bool doCLs, int ntoysM2lnQ, bool doBys, int ntoysBys){
		fAlpha=alpha; fConfidenceLevel=1-fAlpha; _noutcomes=noutcome;
		_doCLs=doCLs; _doBys=doBys; _ntoysM2lnQ=ntoysM2lnQ; _ntoysBys=ntoysBys; 
		if(_clslimit){ _clslimit->SetAlpha(fAlpha); }
		if(_byslimit){_byslimit->SetAlpha(fAlpha); _byslimit->SetNumToys(_ntoysBys);  _noutcomes = _noutcomes;}
		Bands();
	}
	void LimitBands::Bands(){
		clock_t start_time=clock(), cur_time=clock(), funcStart_time=clock();
		if(!_doCLs && !_doBys ) {cout<<" Bands, do nothing"<<endl; return;}
		if(!_cms) {cout<<" Bands, model not set yet"<<endl; return; }
		if(_doCLs && (!_clslimit || !_frequentist) ) {cout<<" Bands, clslimit or frequentist not set yet"<<endl; return; }		
		if(_doBys && !_byslimit) {cout<<" Bands, byslimit not set yet"<<endl; return; }

		_cms->SetSignalScaleFactor(1.);

		double rmeancls=0, rmeanbys=0;

		_vrcls.clear(); _vpcls.clear(); _difrcls.clear();
		_vrbys.clear(); _vpbys.clear(); _difrbys.clear();

		vector<double> vrcls; vrcls.clear();
		vector<double> vpcls; vpcls.clear();
		vector<double> vrbys; vrbys.clear();
		vector<double> vpbys; vpbys.clear();

		if(_debug >= 10){ // print out s,b se, be, correlation, err_pdf_type
			cout<<"\t DoingStatisticalBandsForLimit print out details on signal, background and uncertainties ...."<<endl;
			_cms->Print();
		}

		vector< VIChannel > vvPossibleOutcomes; vvPossibleOutcomes.clear();
		vector<int> vEntries; vEntries.clear();  // index of vvPossibleOutcomes are the same as vEntries
		VIChannel vbkg_tmp; 

		vector< vector< RooAbsData* > > vvPossibleOutcomes_forUnbinnedChannels;
		vvPossibleOutcomes_forUnbinnedChannels.clear();

		if(_debug) cout<<" _noutcomes="<<_noutcomes<<endl;

		double *pars = NULL;
		if(_tossPseudoDataConvention==1){
			pars = _cms->Get_fittedParsInData_b();
			if(pars==NULL) DoAfit(0, _cms->Get_v_data(), _cms->Get_v_pdfs_roodataset(), pars);
		}

		for(int n=0; n<_noutcomes; n++){
			vbkg_tmp.clear(); 
			vbkg_tmp=_cms->GetToyData_H0(pars);
			if(!_cms->hasParametricShape()){
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
			}else{
				vvPossibleOutcomes.push_back(vbkg_tmp);
				vEntries.push_back(1);	
				vector<RooAbsData*> vrds; vrds.clear();
				for(int c=0; c<_cms->Get_vv_pdfs().size(); c++){
					RooDataSet *rds = new RooDataSet(*((RooDataSet*)_cms->Get_v_pdfs_roodataset_toy()[c]));
					vrds.push_back((RooAbsData*)rds);
				}
				vvPossibleOutcomes_forUnbinnedChannels.push_back(vrds);
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

		if(!_cms->hasParametricShape()){
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
		}else{
			for(int i=0; i<npsbl_outcomes; i++){
				iEntries_v[i]=i;
			}
		}

		double nchs=_cms->NumOfChannels();
		int tmpcount;
		if(npsbl_outcomes<50) tmpcount = int(npsbl_outcomes/10);
		else {
			tmpcount = int(npsbl_outcomes/100);
		}

		cout<<"\t Computing expected Upper Limit for "<<_noutcomes<<" sets of outcomes "<<endl;
		double r;

		std::map<double, TTree*> gridCLsb; //r, <clsb, clserr>
		std::map<double, TTree*> gridCLb; //r, <clsb, clserr>
		std::map<double, double> gridQdata; //r, q_data
		if(_doCLs && bM2lnQGridPreComputed){
			ReadM2lnQGridFromFile(sFileM2lnQGrid, gridCLsb, gridCLb, _debug);
			_frequentist->checkFittedParsInData();
		}
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
			_cms->SetData( vvPossibleOutcomes.at(iEntries_v[n]), false ); //set the pseudo data
			if(_cms->hasParametricShape()){
				_cms->SetDataForUnbinned(vvPossibleOutcomes_forUnbinnedChannels[n], false); // set the pseudo data
			}


			int nentries_for_thisR=(int)(vEntries.at(iEntries_v[n]));
			//int nentries_for_thisR=(int)(vEntries.at(iEntries_v[n]));
			if(_doCLs) { 
				if(bOnlyEvalCL_forVR){
					vector<double> v;
					for (int ii=0; ii<_vR_toEval.size(); ii++){
						double rVal = _vR_toEval[ii];
						_cms->SetSignalScaleFactor(rVal);
						if(_frequentist->GetTestStatistics()==1)_frequentist->prepareLogNoverB();
						_frequentist->BuildM2lnQ_data();
						v.push_back( _frequentist->Get_m2lnQ_data() );
					}
					_vvCL_forVR.push_back(v);
					continue;
				}
				//--------------calc r95% with (s,b,n) by throwing pseudo experiments and using the fits to get r95% at CLs=5
				if(!bM2lnQGridPreComputed){
					if(bQuickEstimateInitialLimit) r=_clslimit->LimitOnSignalScaleFactor(_cms, _frequentist, _ntoysM2lnQ);
					else r= _clslimit->LimitOnSignalScaleFactor(_cms, initialRmin, initialRmax, _frequentist, _ntoysM2lnQ, nSteps);

				}else {
					int ii = 0;
					for (std::map<double, TTree *>::iterator itg = gridCLsb.begin(), edg = gridCLsb.end(); itg != edg; ++itg) {
						double rVal = itg->first;
						_cms->SetSignalScaleFactor(rVal);
						if(_frequentist->GetTestStatistics()==1)_frequentist->prepareLogNoverB();
						_frequentist->BuildM2lnQ_data();
						gridQdata[rVal] = _frequentist->Get_m2lnQ_data();
						ii++;
					}
					cout<<" _plotLevel = "<<_plotLevel<<endl;
					TString plotName = "plots_rVScl_toy_"; plotName+=n;
					r = _frequentist->FindLimitFromPreComputedGrid(gridCLsb, gridCLb, gridQdata, fAlpha, _plotLevel<10?"":plotName);
				}
				vrcls.push_back(r);
				vpcls.push_back(p);
				rmeancls+=r*p;
				for(int ntmp=0; ntmp<nentries_for_thisR; ntmp++){
					_difrcls.push_back(r);
				}
				if(_debug) {
					vector<double> vvr, vvcls;
					vvr.clear(); vvcls.clear();
					vvr=_clslimit->GetvTestedScaleFactors();	
					vvcls=_clslimit->GetvTestedCLs();
					cout<<" DoingStatisticalBandsForLimit cls_r95\% for n="<<n<<", printing out the recorded r and cls:"<<endl;
					for(int i=0; i<(int)vvr.size(); i++)
						printf("r %5.2f cls %3.2f\n",vvr[i], vvcls[i]);
				}
				if(_debug)cout<<"n="<<n<<" cls_r= "<<r<<" p= "<<p<<" repeated "<<nentries_for_thisR<<endl;
				if(_debug) { start_time=cur_time; cur_time=clock(); cout << "\t\t\t TIME_in_BAND, doCLs "<< n << " took " <<(cur_time - start_time)/1000000. << " sec\n"; }
				if(_plotLevel>=10){
					TPaveText *pt = SetTPaveText(0.2, 0.7, 0.3, 0.9);
					DrawEvolution2D d2d(_clslimit->GetvTestedScaleFactors(), _clslimit->GetvTestedCLs(), "; r ; CLs", TString::Format("expectation_r_vs_cl_%d", n).Data(), pt);
					d2d.draw();
					d2d.save();
				}
			}
			if(_doBys){
				_cms->SetSignalScaleFactor(1.);
				r=_byslimit->Limit(_cms);
				if(r>0){
					vrbys.push_back(r);
					vpbys.push_back(p);
					rmeanbys+=r*p;
					for(int ntmp=0; ntmp<nentries_for_thisR; ntmp++){
						_difrbys.push_back(r);
					}
				}else{
					_actualOutComesForBys -= (int)(_noutcomes*p);
				}
				if(_debug)cout<<"n="<<n<<" bys_r= "<<r<<" p= "<<p<<" repeated "<<nentries_for_thisR<<"  debug "<<rmeanbys<<endl;
				if(_debug) { start_time=cur_time; cur_time=clock(); cout << "\t\t\t TIME_in_BAND, doBys "<< n << " took " <<(cur_time - start_time)/1000000. << " sec\n"; }
			}
			if(_doCLs && _doBys && _debug ) {
				if(vrcls.back()==0 || vrbys.back()==0) cout<<"Warning:  r=0 "<<endl;
				else {
					cout<<"\n\t CompareCLsBys n="<<n<<" cls= "<<vrcls.back()<<" bys= "<<vrbys.back()
						<<" (c-b)/(c+b)="<<(vrcls.back()-vrbys.back())/(vrcls.back()+vrbys.back())<<"  repeated "<<nentries_for_thisR<<endl;
				}
			}

			if(_debug>=10){
				cout<<"******************* pseudo_outcome ("<<n<<") *****************************"<<endl;
				cout<<"         r = "<<r<<endl;
				if(vvPossibleOutcomes[0].size()!=0){
					cout<<" Counting part: "<<endl;
					for(int ii=0; ii<vvPossibleOutcomes.at(iEntries_v[n]).size(); ii++){
						if(ii%10==0) cout<<endl<<"       ";
						printf(" %5d ", vvPossibleOutcomes.at(iEntries_v[n])[ii]);
					}
					cout<<endl;
				}

				if(_cms->hasParametricShape()){
					cout<<" Unbinned part: "<<endl;
					for(int ii=0; ii<vvPossibleOutcomes_forUnbinnedChannels[n].size(); ii++){
						cout<<" * "<<_cms->Get_v_pdfs_channelname()[ii]<<": "<<endl;
						vvPossibleOutcomes_forUnbinnedChannels[n][ii]->Print("V");	
						for(int ie=0; ie<vvPossibleOutcomes_forUnbinnedChannels[n][ii]->sumEntries(); ie++){
							vvPossibleOutcomes_forUnbinnedChannels[n][ii]->get(ie)->Print("V");	
						}
					}
				}
				cout<<endl;

			}
			fflush(stdout);
		}
		if(bOnlyEvalCL_forVR) return;
		// correct weight of each r for bys 
		if(_doBys) {
			if(_actualOutComesForBys<=0) {cout<<"ERROR: LimitBands bys _actualOutComesForBys="<<_actualOutComesForBys<<endl; exit(0);};
			for(int ii=0; ii<vrbys.size(); ii++){
				vpbys[ii] = (vpbys[ii]*_noutcomes/(double)_actualOutComesForBys);
			}
		}

		if(Entries_v) delete [] Entries_v;
		if(iEntries_v) delete [] iEntries_v;

		if(_doCLs) {
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
				cout<<" ProjectingLimits, printing out for cls the r and cummulative p for all possible outcomes:"<<endl;
				for(int i=0; i<(int)_vrcls.size(); i++)
					printf("r %5.2f p %5.4f\n",_vrcls[i], _vpcls[i]);
			}

			GetBandsByFermiCurveInterpolation(_vrcls,_vpcls, _rcls[1], _rcls[3], _rcls[0], _rcls[4]);
			_rcls[5]=rmeancls; _rcls[2]=GetBandByFermiCurveInterpolation(_vrcls, _vpcls, 0.5);
			if(_debug) {
				cout<<"\t ProjectingLimits projecting cls r at (these are interpolated by Fermi fuction, now should be deprecated)"<<endl; 
				cout<<"    -2 sigma ="<<_rcls[0]<<endl;
				cout<<"    -1 sigma ="<<_rcls[1]<<endl;
				cout<<"    median   ="<<_rcls[2]<<endl;
				cout<<"     1 sigma ="<<_rcls[3]<<endl;
				cout<<"     2 sigma ="<<_rcls[4]<<endl;
				cout<<"     mean    ="<<_rcls[5]<<endl;
			}
		}
		if(_doBys) {
			int nr = vrbys.size(); 
			double *r_v = new double[nr];	
			int *ir_v = new int[nr];	
			for(int i=0; i<nr; i++) r_v[i]=vrbys.at(i);
			Sort(nr, r_v, ir_v, 0); // from small to large
			for(int i=0; i<nr; i++){
				_vrbys.push_back(vrbys.at(ir_v[i]));
				if(i==0) _vpbys.push_back(vpbys.at(ir_v[i]));
				else _vpbys.push_back( _vpbys.at(i-1) + vpbys.at(ir_v[i]) ); 
			}	
			if(r_v) delete [] r_v;
			if(ir_v) delete [] ir_v;

			GetBandsByFermiCurveInterpolation(_vrbys,_vpbys, _rbys[1], _rbys[3], _rbys[0], _rbys[4]);
			_rbys[5]=rmeanbys; _rbys[2]=GetBandByFermiCurveInterpolation(_vrbys, _vpbys, 0.5);
			if(_debug) {
				cout<<"\t ProjectingLimits projecting bys r at (these are interpolated by Fermi fuction, now should be deprecated)"<<endl; 
				cout<<"    -2 sigma ="<<_rbys[0]<<endl;
				cout<<"    -1 sigma ="<<_rbys[1]<<endl;
				cout<<"    median   ="<<_rbys[2]<<endl;
				cout<<"     1 sigma ="<<_rbys[3]<<endl;
				cout<<"     2 sigma ="<<_rbys[4]<<endl;
				cout<<"     mean    ="<<_rbys[5]<<endl;
			}
		}

		if( _doCLs && _doBys ) {
			if(_rcls[5]==0 || _rbys[5]==0) cout<<"Warning:  rmean=0 "<<endl;
			else {
				cout<<"\n\t CompareCLsBys rmean cls= "<<_rcls[5]<<" bys= "<<_rbys[5]
					<<" (c-b)/(c+b)="<<(_rcls[5]-_rbys[5])/(_rcls[5]+_rbys[5])<<endl;
			}
		}
		if(_debug) { start_time=cur_time; cur_time=clock(); cout << "\t\t\t TIME_in_BAND, final printing took " <<(cur_time - start_time)/1000. << " milisec\n"; }
		if(_debug) { cur_time=clock(); cout << "\t\t\t TIME_in_BAND finishing "<<_noutcomes<<" outcomes: " << (cur_time - funcStart_time)/1000000./60. << " minutes\n"; }



		for(int i=0; i<vvPossibleOutcomes_forUnbinnedChannels.size(); i++){
			for(int j=0; j<vvPossibleOutcomes_forUnbinnedChannels[i].size(); j++){
				//		delete (vvPossibleOutcomes_forUnbinnedChannels[i][j]);
			}
		}
	}	
};
