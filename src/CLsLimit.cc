/*
 * =====================================================================================
 * CMS 
 * =====================================================================================
 */
#include <iostream>
#include <iomanip>
#include <math.h>
#include <ctime> // upto second
#include <time.h> // upto micro second

#include "CLsLimit.h"
#include "CRandom.h"
#include "PdfRandom.h"
#include "Utilities.h"
#include "BayesianBase.h"

using std::cout;
using std::endl;
using std::min;
namespace lands{

	// Class CLsBase
	CLsBase::CLsBase(){
		Q_b = 0;
		Q_sb = 0;
		iq_sb=0;
		iq_b=0;
		_nsig=0; _nbkg=0; _ndat=0;
		_debug=0;
		_rdm=0;
	}

	CLsBase::~CLsBase(){
		if(Q_b) delete [] Q_b;
		if(Q_sb) delete [] Q_sb;
		if(iq_b)delete [] iq_b;
		if(iq_sb)delete [] iq_sb;
		_rdm=0;
	}
	bool CLsBase::BuildM2lnQ(CountingModel *cms, int nexps, int sbANDb_bOnly_sbOnly){
		_model=cms;
		BuildM2lnQ(nexps, sbANDb_bOnly_sbOnly);
	}
	bool CLsBase::BuildM2lnQ(int nexps, int sbANDb_bOnly_sbOnly){  // 0 for sbANDb, 1 for bOnly, 2 for sbOnly
		if(!_model) { 
			cout<<"No model constructed....exit"<<endl;
			exit(0);
		}
		if(! (_model->Check()) ){
			cout<<"Model is not correctly constructed, exit"<<endl;
			_model->Print();
			exit(0);
		}
		_rdm=_model->GetRdm();

		clock_t start_time=clock(), cur_time=clock();

		if( _debug >= 100 )_model->Print();

		//------if input is null, then do nothing	
		_nexps = nexps;
		_nchannels = _model->NumOfChannels();

		_nsig=0; _nbkg=0; _ndat=0;
		vector<double> vs, vb, vd; 
		vs.clear(); vb.clear(); vd.clear();
		for(int i=0; i<_nchannels; i++){
			double totbkg = 0;
			for(int isamp = 1; isamp<(_model->Get_vv_exp_sigbkgs())[i].size(); isamp++){
				totbkg+=(_model->Get_vv_exp_sigbkgs())[i][isamp];
			}
			vs.push_back((_model->Get_vv_exp_sigbkgs())[i][0]);
			vb.push_back(totbkg);
			vd.push_back((_model->Get_v_data())[i]);
			_nsig += vs[i];
			_nbkg += vb[i];
			_ndat += vd[i];
		}

		if(Q_sb) delete [] Q_sb;
		if(Q_b) delete [] Q_b;
		if(iq_sb) delete [] iq_sb;
		if(iq_b) delete [] iq_b;
		Q_sb=new double[_nexps];
		Q_b=new double[_nexps];
		iq_sb = new int[_nexps];	
		iq_b = new int[_nexps];	

		Q_b_exp = 0; //_nbkg*log((_nsig+_nbkg)/_nbkg);
		Q_b_data     = 0;


		double *n=new double[_nchannels];
		double *noverb=new double[_nchannels];
		double *lognoverb=new double[_nchannels];
		for(int i=0; i<_nchannels; i++){	
			// skip a channle in which nsig==0 || ntotbkg==0
			if(vs[i] > 0 && vb[i] > 0) {
				n[i]=vs[i]+vb[i];
				noverb[i]=n[i]/vb[i];
				lognoverb[i]=log(noverb[i]);
				Q_b_exp+=(vb[i]*lognoverb[i]);
				Q_b_data    +=(vd[i]*lognoverb[i]);
			}else{
				lognoverb[i]=0;
			}

			if(_debug>=10)cout<<" \t channel "<<i<<" s="<<vs[i]<<" b="<<vb[i]<<" d="<<vd[i]<<endl;	

		}
		int nsbi, nbi;
		int tenth = _nexps/10;
		int ntemp = _nexps*_nchannels;
		if(_debug >=10 ) cout<<"ntemp="<<ntemp<<endl;
		if( ntemp>=10000000 ) {
			cout<<"\t gonna generate "<<ntemp*2<<" poisson numbers "<<endl;
		}
		for(int i=0; i<_nexps; i++){
			if( ntemp>=10000000 ) {
				if( (i+1)%tenth == 0 ){
					printf("... Building -2lnQ,  %4.1f \%\n", i/(double)_nexps*100);
					fflush(stdout);
				}
			}
			Q_sb[i]=0;Q_b[i]=0;	
			vector< vector<double> > vv = _model->FluctuatedNumbers();
			for(int ch=0; ch<_nchannels; ch++){	
				double totbkg = 0; 
				for(int isamp=1; isamp<vv[ch].size(); isamp++){
					totbkg+=vv[ch][isamp];
				}
				if(vv[ch][0] > 0 && totbkg > 0 ) {
					if( sbANDb_bOnly_sbOnly != 1 ){
						nsbi = _rdm->Poisson( vv[ch][0] + totbkg );	
						Q_sb[i] += (nsbi*lognoverb[ch]) ;
					}
					if( sbANDb_bOnly_sbOnly != 2 ){
						nbi  = _rdm->Poisson(totbkg);
						Q_b[i]  += (nbi *lognoverb[ch]) ;
					}
				}
			}
		}

		if(_debug) { start_time=cur_time; cur_time=clock(); cout << "\t\t\t TIME in RunMCExps run_"<<_nexps<<"_pseudo exps: " << (cur_time - start_time)/1000. << " millisec\n"; }

		//ProcessM2lnQ();
		if( sbANDb_bOnly_sbOnly !=2 )Sort(_nexps, Q_b, iq_b, 0); // rank from small to large
		if( sbANDb_bOnly_sbOnly !=1 )Sort(_nexps, Q_sb, iq_sb, 0);
		if (sbANDb_bOnly_sbOnly==0) 
			if( ( Q_b_data < Q_b[iq_b[0]] || Q_b_data > Q_sb[iq_sb[_nexps-1]] ) ){ 
				cout<<"\t probability of this -2lnQ_data is very very small, it's out of "<<_nexps<<" exps, you need more toys"<<endl;
				cout<<"\t -2lnQ_data =   "<<-2*Q_b_data+2*_nsig<<endl;
				cout<<"\t -2lnQ_b    = [ "<<-2*Q_b[iq_b[_nexps-1]]+2*_nsig<<" , "<<-2*Q_b[iq_b[0]]+2*_nsig<<" ]"<<endl;
				cout<<"\t -2lnQ_sb   = [ "<<-2*Q_sb[iq_sb[_nexps-1]]+2*_nsig<<" , "<<-2*Q_sb[iq_sb[0]]+2*_nsig<<" ]"<<endl;
			}
		if(_debug>=1 && sbANDb_bOnly_sbOnly==0 ) {
			cout<<"\t -2lnQ_data =   "<<-2*Q_b_data+2*_nsig<<endl;
			cout<<"\t -2lnQ_b    = [ "<<-2*Q_b[iq_b[_nexps-1]]+2*_nsig<<" , "<<-2*Q_b[iq_b[0]]+2*_nsig<<" ]"<<endl;
			cout<<"\t -2lnQ_sb   = [ "<<-2*Q_sb[iq_sb[_nexps-1]]+2*_nsig<<" , "<<-2*Q_sb[iq_sb[0]]+2*_nsig<<" ]"<<endl;
		}
		if(_debug>=100 && sbANDb_bOnly_sbOnly==0 ){
			cout<<"\t\t CHECKING order ----- "<<endl;
			for(int i=0; i<_nexps; i++){
				cout<<"\t "<<i<<"     "<<-2*Q_b[iq_b[i]]+2*_nsig<<"  "<<-2*Q_sb[iq_sb[i]]+2*_nsig<<endl;
			}
		}

		delete [] n; delete [] noverb; delete [] lognoverb;
		return true;
	}

	void CLsBase::ProcessM2lnQ(){
		Sort(_nexps, Q_b, iq_b, 0); // rank from small to large
		Sort(_nexps, Q_sb, iq_sb, 0);
		if( Q_b_data < Q_b[iq_b[0]] || Q_b_data > Q_sb[iq_sb[_nexps-1]] ){ 
			cout<<"\t probability of this -2lnQ_data is very very small, it's out of "<<_nexps<<" exps, you need more toys"<<endl;
			cout<<"\t -2lnQ_data =   "<<-2*Q_b_data+2*_nsig<<endl;
			cout<<"\t -2lnQ_b    = [ "<<-2*Q_b[iq_b[_nexps-1]]+2*_nsig<<" , "<<-2*Q_b[iq_b[0]]+2*_nsig<<" ]"<<endl;
			cout<<"\t -2lnQ_sb   = [ "<<-2*Q_sb[iq_sb[_nexps-1]]+2*_nsig<<" , "<<-2*Q_sb[iq_sb[0]]+2*_nsig<<" ]"<<endl;
		}
		if(_debug>=1) {
			cout<<"\t -2lnQ_data =   "<<-2*Q_b_data+2*_nsig<<endl;
			cout<<"\t -2lnQ_b    = [ "<<-2*Q_b[iq_b[_nexps-1]]+2*_nsig<<" , "<<-2*Q_b[iq_b[0]]+2*_nsig<<" ]"<<endl;
			cout<<"\t -2lnQ_sb   = [ "<<-2*Q_sb[iq_sb[_nexps-1]]+2*_nsig<<" , "<<-2*Q_sb[iq_sb[0]]+2*_nsig<<" ]"<<endl;
		}
		if(_debug>=100){
			cout<<"\t\t CHECKING order ----- "<<endl;
			for(int i=0; i<_nexps; i++){
				cout<<"\t "<<i<<"     "<<-2*Q_b[iq_b[i]]+2*_nsig<<"  "<<-2*Q_sb[iq_sb[i]]+2*_nsig<<endl;
			}
		}
	}
	void CLsBase::SetRdm(CRandom *rdm){
		_rdm = rdm;
	}
	vector<double> CLsBase::Get_m2logQ_sb(){
		vector<double> tmp;
		for(int i=0; i<_nexps; i++){
			tmp.push_back(-2*(Q_sb[i]-_nsig));
		}
		return tmp;
	} 
	vector<double> CLsBase::Get_m2logQ_b(){
		vector<double> tmp;
		for(int i=0; i<_nexps; i++){
			tmp.push_back(-2*(Q_b[i]-_nsig));
		}
		return tmp;
	} 
	double CLsBase::CLsb(){
		double ret =0;// 1./(double)_nexps;
		double tmp = Q_b_data;
		for(int i=0; i<_nexps; i++){
			if(Q_sb[iq_sb[i]]<=tmp)
				ret = (i+1)/(double)_nexps;	
		}		

		if(_debug>=10){
			cout<<"CLsBase::CLsb  CLsb()="<<ret<<" and total exps="<<_nexps<<endl;
			if(ret*_nexps <= 20) cout<<"CLsBase::CLsb  CLsb*nexps="<<ret*_nexps<<", statistic may not enough"<<endl;
		}
		if(ret == 0){
			cout<<"CLsBase::CLsb CLsb=0, it means number of pseudo experiments is not enough"<<endl;
			cout<<"              Currently, we put CLsb=1./"<<_nexps<<endl;
			ret = 1./(double)_nexps;
		}
		return ret;
	}
	double CLsBase::CLs(){
		double clb=CLb();
		double clsb=CLsb();
		if(clb==0){cout<<"CLsBase::CLs  Warning clb_b==0 !!!!"<<endl; return 1;}
		if(_debug>=10) cout<<"CLsBase::CLs  CLs=CLsb/CLb="<<clsb/clb<<endl;
		return clsb/clb;
	}
	double CLsBase::CLb(){
		return CLb(Q_b_data);
	}
	double CLsBase::PValue(double lnq){
		double ret=0;
		bool hasQ_gt_lnq = false;
		for(int i=0; i<_nexps; i++){ 
			if(Q_b[iq_b[i]] >= lnq)  {
				ret = i/(double)_nexps;	
				hasQ_gt_lnq=true;
				break;
			}
		}		
		if(hasQ_gt_lnq==false) {
			ret= 1-1./(double)_nexps;
			/*
			   cout<<"********WARNING********"<<endl;
			   cout<<" Toys for b-only hypothesis  may be not enough to evaluate the true significance, "<<endl;
			   cout<<" Q_b[0]="<<Q_b[iq_b[0]]<<" Q_b["<<_nexps<<"]="<<Q_b[iq_b[_nexps-1]]
			   <<", and tested Q="<<lnq<<endl;	
			   cout<<" we set PValue to be 1./_nexps = "<<1-ret<<endl;
			 */
		}
		return 1-ret;
	}
	void CLsBase::CheckFractionAtHighEnd(vector<double> vlogQ, vector<double> vlogQ_prob){
		//cout<<endl<<"*********Start Calc mean value of significance ......."<<endl;
		// vlogQ has been sorted from small to larger ...,  vlogQ_prob for accumulative probability
		double fractionGTmaxlnQb = 0;
		for (int i=0; i<vlogQ.size(); i++) {
			if(vlogQ[i] > Q_b[iq_b[_nexps-1]]) {
				if(i>0) fractionGTmaxlnQb = 1 - vlogQ_prob[i-1];	
				else fractionGTmaxlnQb = 1;
				break;
			}
		}
		cout<<" Fraction of logQ from S+B hypothesis larger than maximum logQ_b  = "<<fractionGTmaxlnQb<<endl;

	}
	double CLsBase::CLb(double lnq){
		double ret =0;// 1./(double)_nexps;
		double tmp = lnq;
		for(int i=0; i<_nexps; i++){ 
			if(Q_b[iq_b[i]]<=tmp)
				ret = (i+1)/(double)_nexps;	
		}		
		if(_debug>=10){
			cout<<"CLsBase::CLb  CLb()="<<ret<<" and total exps="<<_nexps<<endl;
			int step = _nexps/20;
			cout<<"*** print out Q_b and accumulative probability, size="<<_nexps<<" step="<<step<<endl;
			int i=0;
			for(; i<_nexps; i+=(step+1)) {
				printf("%10d",i);
				cout<<"\t Q_b,p= "<<Q_b[iq_b[i]]<<" "<<double((i+1)/(double)_nexps)<<endl;
			}
			if(i!=_nexps || (i==_nexps && step!=1)) {
				printf("%10d",_nexps);
				cout<<"\t lnQ,p= "<<Q_b[iq_b[_nexps-1]]<<" 1"<<endl;
			}
		}
		if(ret*_nexps <= 20) cout<<"CLsBase::CLb  CLb*nexps="<<ret*_nexps<<", statistic may not enough"<<endl;
		if( (1-ret)*_nexps <= 20) cout<<"CLsBase::CLb  (1-CLb)*nexps="<<(1-ret)*_nexps<<", statistic may not enough"<<endl;
		if(ret == 0 || ret==1){
			cout<<"CLsBase::CLb CLb="<<ret<<", it means number of pseudo experiments is not enough"<<endl;
			if(ret==0) {
				cout<<"              Currently, we put CLb=1./"<<_nexps<<endl;
				ret = 1./(double)_nexps;
			}
			if(ret==1) {
				cout<<"              Currently, we put CLb= 1 - 1./"<<_nexps<<endl;
				ret = 1 - 1./(double)_nexps;
			}
		}
		return ret;
	}
	double CLsBase::CLsb_b(){
		double ret = 1./(double)_nexps;
		double tmp = Q_b_exp;
		for(int i=0; i<_nexps; i++){
			if(Q_sb[iq_sb[i]]<=tmp)
				ret = (i+1)/(double)_nexps;	
		}		
		return ret;
	}
	double CLsBase::CLs_b(){
		double clb_b=CLb_b();
		double clsb_b=CLsb_b();
		if(clb_b==0){cout<<"Warning clb_b==0 !!!!"<<endl;return 1;}
		return clsb_b/clb_b;
	}
	double CLsBase::CLb_b(){
		double ret = 1./(double)_nexps;
		double tmp = Q_b_exp;
		for(int i=0; i<_nexps; i++){
			if(Q_b[iq_b[i]]<=tmp)
				ret = (i+1)/(double)_nexps;	
		}		
		return ret;
	}
	void CLsBase::SetLogQ_b(vector<double> vlnQ_b){
		_nexps=vlnQ_b.size();
		Q_b=new double[_nexps];
		iq_b = new int[_nexps];	
		for(int i=0; i<_nexps; i++) Q_b[i]=vlnQ_b[i];
		Sort(_nexps, Q_b, iq_b, 0);
	}
	void CLsBase::SetLogQ_sb(vector<double> vlnQ_sb){
		_nexps=vlnQ_sb.size();
		Q_sb=new double[_nexps];
		iq_sb = new int[_nexps];	
		for(int i=0; i<_nexps; i++) Q_sb[i]=vlnQ_sb[i];
		Sort(_nexps, Q_sb, iq_sb, 0);
	}
	void CLsBase::SetLogQ_data(double lnQ_data){Q_b_data=lnQ_data;}

	double CLsBase::Get_m2lnQ_data(){return -2*(Q_b_data-_nsig);}
	void CLsBase::SetDebug(int debug){_debug=debug;}
	CRandom* CLsBase::GetRdm(){return _rdm;}

	void CLsBase::tmpFun0(vector<double> & vlogQ, vector<double>& vlogQ_prob){
		SortAndCumulative(Q_sb, _nexps, vlogQ, vlogQ_prob, 0);// sort it by  increased order 
	}
	double CLsBase::SignificanceComputation(int ntoys_for_sb, int ntoys_for_b){
		vector<double> vsignificance, vsignificance_cp;
		return SignificanceComputation(ntoys_for_sb, ntoys_for_b, vsignificance, vsignificance_cp);	
	}
	double CLsBase::SignificanceComputation(int ntoys_for_sb, int ntoys_for_b, vector<double>& vsignificance, vector<double> & vsignificance_cp){
		clock_t start_time, cur_time, funcStart_time;
		start_time=clock(); cur_time=clock(); funcStart_time=clock();

		if(_debug>=10)cout<<" previous _nexps = "<<_nexps<<endl;

		vector<double> vlogQ_sb, vlogQ_sb_prob;
		vlogQ_sb.clear(); vlogQ_sb_prob.clear(); vsignificance_cp.clear(); vsignificance.clear();

		if(ntoys_for_sb<=0) {
			if(_debug>=10)cout<<" using the old _nexps for sb = "<<_nexps<<endl;
			SortAndCumulative(Q_sb, _nexps, vlogQ_sb, vlogQ_sb_prob, 0);// sort it by  increased order 
		}else{
			if(_debug>=10)cout<<" producing new _nexps for sb = "<<ntoys_for_sb<<endl;
			BuildM2lnQ(ntoys_for_sb, 2);
			if(_debug>=10)cout<<" end new _nexps for sb = "<<ntoys_for_sb<<endl;
			SortAndCumulative(Q_sb, _nexps, vlogQ_sb, vlogQ_sb_prob, 0);// sort it by  increased order 
			if(_debug>=10)cout<<" sort _nexps for sb = "<<ntoys_for_sb<<endl;
			if(_debug){
				start_time=cur_time; cur_time=clock(); cout << "\t TIME_in_BuildM2logQsb: "<< ntoys_for_sb <<"toys, "<< (cur_time - start_time)/1000000. << " sec\n";
				fflush(stdout);
			}
		}	

		if(_debug>=10)cout<<" producing new _nexps for bonly = "<<ntoys_for_b<<endl;
		BuildM2lnQ(ntoys_for_b, 1);
		if(_debug){
			start_time=cur_time; cur_time=clock(); cout << "\t TIME_in_BuildM2logQb: "<< ntoys_for_b <<"toys, "<< (cur_time - start_time)/1000000./60. << " minutes\n";
			fflush(stdout);
		}
		CheckFractionAtHighEnd(vlogQ_sb, vlogQ_sb_prob);

		int previousStopPoint=0;
		for(int isb=0; isb<vlogQ_sb.size(); isb++) {
			if(_debug>=100)cout<<"previousStopPoint = "<<previousStopPoint<<endl;
			double lnq=vlogQ_sb[isb];
			double ret=0;
			bool hasQ_gt_lnq = false;
			for(int i=previousStopPoint; i<_nexps; i++){ 
				if(Q_b[iq_b[i]] >= lnq)  {
					ret = i/(double)_nexps;	
					hasQ_gt_lnq=true;
					previousStopPoint = i;
					break;
				}
			}		
			if(hasQ_gt_lnq==false) {
				ret= 1-1./(double)_nexps;
			}

			double tmp = 1-ret;
			if(tmp>0.5) vsignificance.push_back(0);
			else vsignificance.push_back( Significance(tmp) );		
		}

		if(_debug){
			start_time=cur_time; cur_time=clock(); cout << "\t TIME_in_sortsigs: " << (cur_time - start_time)/1000000. << " secs\n"; fflush(stdout);
		}
		vsignificance_cp=vlogQ_sb_prob;
		double significance_mean = GetMeanOfSortedXwithProb(vsignificance, vlogQ_sb_prob);
		return significance_mean;
	}
	double CLsBase::SignificanceForData(int ntoys_for_b){
		if(ntoys_for_b > 0)  {
			BuildM2lnQ(ntoys_for_b, 1);
		}
		double pvalue=PValue(Q_b_data);
		double significance = Significance(pvalue);
		return significance;
	}

	double CLsBase::lnQsb_sigma(int sigma){
		double ret;	
		vector<double> vlogQ, vlogQ_prob;
		SortAndCumulative(Q_sb, _nexps, vlogQ, vlogQ_prob, 0);// sort it by  increased order 
		if(sigma<=2 && sigma>=-2){
			double r[5];	
			if(sigma!=0){
				GetBandsByLinearInterpolation(vlogQ,vlogQ_prob, r[1], r[3], r[0], r[4] );
				ret=r[sigma+2];
			}
			if(sigma==0)
				ret=GetBandByLinearInterpolation(vlogQ, vlogQ_prob, 0.5);
		}
		else {
			ret=GetMeanOfSortedXwithProb(vlogQ, vlogQ_prob);
		}

		if( _debug>=10 ){
			int step = vlogQ.size()/20;
			cout<<"CLsBase: lnQsb_sigma"<<endl;
			cout<<"*** print out lnQ and accumulative probability, size="<<vlogQ.size()<<" step="<<step<<endl;
			int i=0;
			for(; i<vlogQ.size(); i+=(step+1)) {
				printf("%10d",i);
				cout<<"\t lnQ,p= "<<vlogQ[i]<<" "<<vlogQ_prob[i]<<endl;
			}
			if(i!=vlogQ.size() || (i==vlogQ.size() && step!=1)) {
				printf("%10d",vlogQ.size());
				cout<<"\t lnQ,p= "<<vlogQ.back()<<" "<<vlogQ_prob.back()<<endl;
			}
		}
		return ret;
	}



	// Class CLsLimit	
	double CLsLimit::LimitOnSignalScaleFactor(CountingModel *cms,
			double minRtoScan, double maxRtoScan,
			CLsBase *frequentist, int nexps, int nstep ){
		_frequentist=frequentist; _nexps=nexps; 

		clock_t start_time, cur_time;
		start_time=clock(); cur_time=clock();

		double epsilon = _clstolerance;
		//	int nsigma=5;

		_vR.clear(); _vCLs.clear();

		if(_debug) cout<<"CLsLimit::LimitOnSignalScaleFactor  looking for C.L. 95% Limit on the ratio ----"<<endl;


		double r0,r1;
		double cl0, cl1;

		if(minRtoScan!=maxRtoScan) {

			if(minRtoScan>=maxRtoScan || minRtoScan <=0 ) {
				cout<<"Error in LimitOnSignalScaleFactor: (minRtoScan="<<minRtoScan<<") >= (maxRtoScan="<<maxRtoScan<<", exit"<<endl;
				cout<<"please make sure maxRtoScan > minRtoScan > 0 "<<endl;
				exit(0);
			}
			if(nstep<1)  {cout<<" steps in autoscan should not less than 1, exit"<<endl; exit(0);}
			if(_debug) cout<<"\t First auto scaning R from  "<<minRtoScan<<" to "<<maxRtoScan<<" in "<<nstep<<" steps"<<endl;

			for(double rmid=minRtoScan; rmid<=maxRtoScan; rmid+=(maxRtoScan-minRtoScan)/(double)nstep){
				cms->SetSignalScaleFactor(rmid);
				frequentist->BuildM2lnQ(cms, nexps); 
				cl0=frequentist->CLs(); //_nsigma(nsigma);
				_vR.push_back(rmid);_vCLs.push_back(cl0);
				if(_debug)cout<<"TESTED r="<<rmid<<"  CLs="<<cl0<<endl;
			}
			// --- -get the _r95 @ alpha=0.05 by linear interpolation
			double x1=0, y1=0, x2=0, y2=0;
			int nsize=_vCLs.size();
			double *dCLs=new double[nsize];
			for(int i=0; i<nsize; i++){
				dCLs[i]=_vCLs[i];
			}
			int *ix=new int[nsize];
			Sort(nsize, dCLs, ix, 0);
			for(int i=0; i<nsize; i++){
				if(dCLs[ix[i]]==_alpha) {_r95=_vR[ix[i]]; return _r95;}
				if(dCLs[ix[i]]<_alpha) { 
					x1=_vR[ix[i]]; 
					y1=dCLs[ix[i]];
				}
				if(dCLs[ix[i]]>_alpha && x2==0){
					x2=_vR[ix[i]]; 
					y2=dCLs[ix[i]];
				}
			}
			if(!x1 || !x2 || !y1 || !y2){
				cout<<"Warning: Your initial estimated R range is not suitable, all CLs values I got are in one side of "<<_alpha<<endl;
				x1 = _vR[0]; x2=_vR.back();
				y1 = _vCLs[0]; y2=_vCLs.back();
			}
			if(ix) delete []ix;
			if(dCLs) delete [] dCLs;
			r0=x1; r1=x2; cl0=y1; cl1=y2;
		} else {
			cms->SetSignalScaleFactor(1.);
			//		if(_debug>=10)
			//			r0=R_CL(cms, 1-_alpha, 1.e-30, 1.e-6, 1000., _debug );
			//		else 
			//			r0=R_CL(cms, 1-_alpha, 1.e-30, 1.e-6, 1000., 0);
			//		if(r0>100) r0=10; // FIXME will use the new Bayesian estimation
			BayesianBase bys(cms, 0.05, 1.e-2);
			bys.SetNumToys(100);
			bys.SetDebug(_debug);
			r0= bys.Limit();
			if(_debug){
				start_time=cur_time; cur_time=clock(); cout << "\t TIME_in_BayesianLimitEsitmate: "<< (cur_time - start_time) << " microsec\n";
			}

			if(_debug)cout<<"CLsLimit::LimitOnSignalScaleFactor: start with estimate from Bayesian technique, r95%="<<r0<<endl;

			cms->SetSignalScaleFactor(r0);
			frequentist->BuildM2lnQ(cms, nexps); 
			cl0=frequentist->CLs();//_nsigma(nsigma);
			_vR.push_back(r0);_vCLs.push_back(cl0);
			if(_debug)cout<<"Estimated_initial r="<<r0<<"  CLs="<<cl0<<endl;
			if(fabs(cl0-_alpha)<=epsilon) {
				_r95=r0;
				if(_debug)cout<<"Converge at CLs="<<_alpha<<"+/-"<<epsilon<<" by "<<"1 iteration"<<endl;
				return r0;
			}

			//	r1=r0*0.90; //----------usually, CLs-limit is more aggresive than Bayesian's, about 10% smaller.
			if(cl0>_alpha) r1=r0*1.10;	
			else r1=r0*(cl0/_alpha);
			//else r1=r0*0.90;

			cms->SetSignalScaleFactor(r1);
			frequentist->BuildM2lnQ(cms, nexps); 
			cl1=frequentist->CLs();//_nsigma(nsigma);
			_vR.push_back(r1);_vCLs.push_back(cl1);
			if(_debug)cout<<"Estimated_r="<<r1<<"  CLs="<<cl1<<endl;
			if(fabs(cl1-_alpha)<=epsilon) {
				_r95=r1;
				if(_debug)cout<<"Converge at CLs="<<_alpha<<"+/-"<<epsilon<<" by "<<"2 iterations"<<endl;
				return r1;
			}
		}

		bool foundit=false;
		double rmid=0;
		int nmaxrepeat=0;	
		while(!foundit && nmaxrepeat<=30 ){ 
			// --- -  --
			//  using linear interpolation to do converge will be quicker
			// ---------
			rmid=LogLinearInterpolation(r0,cl0,r1,cl1,_alpha);
			if(rmid<0) rmid=-rmid;
			if(_debug >= 10 )cout<<" r0="<<r0<<" cl0="<<cl0<<" r1="<<r1<<" cl1="<<cl1<<" rmid="<<rmid<<endl;
			cms->SetSignalScaleFactor(rmid);
			frequentist->BuildM2lnQ(cms, nexps); 
			double clmid=frequentist->CLs();//_nsigma(nsigma);
			_vR.push_back(rmid);_vCLs.push_back(clmid);
			if(_debug)cout<<"TESTED r="<<rmid<<"  CLs="<<clmid<<endl;
			if(fabs(clmid-_alpha)<epsilon){
				foundit=true; //alpha=0.05 C.L. 95% 		
				_r95=rmid;
				if(_debug)cout<<"Converge at CLs="<<_alpha<<"+/-"<<epsilon<<" by "<<_vR.size()<<" iterations"<<endl;
				return rmid;
			}
			else {
				double x[3]={r0,rmid,r1}; double y[3]={cl0,clmid,cl1};		
				int iy[3];
				Sort(3,y,iy,0);
				if(_alpha<y[iy[1]]){ //---------kick out a number among r0, r1, rmid
					cl0=y[iy[0]]; cl1=y[iy[1]];
					r0=x[iy[0]]; r1=x[iy[1]];
				}
				else if(_alpha>y[iy[1]]){
					cl0=y[iy[1]]; cl1=y[iy[2]];
					r0=x[iy[1]]; r1=x[iy[2]];
				}
				else {
					cout<<"Warning: CANNOT converge: r0="<<r0<<" cl0="<<cl0<<" r1="<<r1<<" cl1="<<cl1<<" rmid="<<rmid<<" clmid="<<clmid<<endl;
				}
			}
			nmaxrepeat++;
			//------------------------????????????? is CLs mono-increased/decreased as r ?  has to investigate it 
			if(!foundit){
				if(rmid==_vR[_vR.size()-2]) {
					foundit=true; 
					if(_debug)cout<<"We get rmid="<<rmid<<" twice, and decide to stop here...."<<endl;
				}
			}
			if(foundit && _vR.size()>1){	
				bool hasCLsGT05=false;
				bool hasCLsLT05=false;
				int nmax_tmp=0; 

				// --- find the rmid  with smallest distance btw its cls and 0.05
				double tmp_min_dist = 9999;
				//int tmp_min_index = 0;
				for(int i=0; i<_vCLs.size(); i++ ){
					if( fabs( _vCLs.at(i)-_alpha ) < tmp_min_dist) {
						tmp_min_dist=fabs(_vCLs.at(i)-_alpha);
						rmid=_vR.at(i);
					}
				}	
				if(_debug)cout<<"\t temply we found r="<<rmid<<" with CLs - alpha = "<<tmp_min_dist<<endl;

				while( (!hasCLsLT05 || !hasCLsGT05) && nmax_tmp<30 ) {
					for(int icls= 0; icls<_vR.size(); icls++){
						if(_vCLs[icls]<=_alpha) hasCLsLT05=true;
						if(_vCLs[icls]>=_alpha) hasCLsGT05=true;
					}	
					if(!hasCLsLT05 && _debug) cout<<"\t we don't have CLs < "<<_alpha<<endl;
					if(!hasCLsGT05 && _debug) cout<<"\t we don't have CLs > "<<_alpha<<endl;
					if(!hasCLsLT05) {
						rmid= rmid*1.05; //FIXME this number should more smart 
						cms->SetSignalScaleFactor(rmid);
						frequentist->BuildM2lnQ(cms, nexps); 
						double clmid=frequentist->CLs();//_nsigma(nsigma);
						_vR.push_back(rmid);_vCLs.push_back(clmid);
						if(_debug)cout<<"TESTED r="<<rmid<<"  CLs="<<clmid<<endl;
					}
					if(!hasCLsGT05) {
						rmid= rmid*0.95; // FIXME this number should be more smart
						cms->SetSignalScaleFactor(rmid);
						frequentist->BuildM2lnQ(cms, nexps); 
						double clmid=frequentist->CLs();//_nsigma(nsigma);
						_vR.push_back(rmid);_vCLs.push_back(clmid);
						if(_debug)cout<<"TESTED r="<<rmid<<"  CLs="<<clmid<<endl;
					}
					if(!hasCLsGT05 || !hasCLsLT05 )nmax_tmp++;
				}//while
				if( nmax_tmp > 0 && _debug )cout<<" \t to get both values GT_and_LT_"<<_alpha<<", we tried "<<nmax_tmp<<" more times"<<endl;
			}//foundit
		}
		int nsize=_vR.size();

		if(!foundit){cout<<"Warning: Not converge to "<<_alpha<<"+/-"<<epsilon<<" with "<<nsize<<" iterations"<<endl;}
		else { if(_debug)cout<<"Converge at alpha="<<_alpha<<"+/-"<<epsilon<<" by "<<nsize<<" iterations"<<endl;}	

		if(nsize==1) {_r95=rmid; return rmid;}

		// --- -get the _r95 @ alpha=0.05 by linear interpolation
		double x1=0, y1=0, x2=0, y2=0;

		double *dCLs=new double[nsize];
		for(int i=0; i<nsize; i++){
			dCLs[i]=_vCLs[i];
		}

		int *ix=new int[nsize];
		Sort(nsize, dCLs, ix, 0);
		for(int i=0; i<nsize; i++){
			if(dCLs[ix[i]]==_alpha) {_r95=_vR[ix[i]]; return _r95;}
			if(dCLs[ix[i]]<_alpha) { 
				x1=_vR[ix[i]]; 
				y1=dCLs[ix[i]];
			}
			if(dCLs[ix[i]]>_alpha && x2==0){
				x2=_vR[ix[i]]; 
				y2=dCLs[ix[i]];
			}
		}
		if(!x1 || !x2 || !y1 || !y2)
		{
			cout<<"Warning.. failed to find both points which are supposed to close to "<<_alpha<<".  r1="
				<<x1<<" CLs1="<<y1<<" r2="<<x2<<" CLs2="<<y2<<endl;
			cout<<"listing all tested r and its corresponding CLs :"<<endl;
			for(int i=0; i<_vR.size(); i++) cout<<"  "<<i<<"  r="<<_vR[i]<<" CLs="<<_vCLs[i]<<endl;
			cout<<"Will use the R with CLs closest to "<<_alpha<<endl;
		}
		// --- find the rmid  with smallest distance btw its cls and 0.05
		double tmp_min_dist = 9999;
		int tmp_min_index = 0;
		for(int i=0; i<_vCLs.size(); i++ ){
			if( fabs( _vCLs.at(i)-_alpha ) < tmp_min_dist) {
				tmp_min_dist=fabs(_vCLs.at(i)-_alpha);
				rmid=_vR.at(i);
				tmp_min_index=i;
			}
		}	
		if(_debug) cout<<"\t closest_to_alpha r="<<rmid<<" cls="<<_vCLs.at(tmp_min_index)<<endl;

		if(x1==x2 && _debug ) cout<<"\t Warning: r1=r2, we are using the same r to do interpolation, just because we get it twice and one for CLs>alpha, and the other for CLs<alpha"<<endl;

		_r95=LogLinearInterpolation(x1,y1,x2,y2,_alpha);
		if(_debug) cout<<"CLsLimit::LimitOnSignalScaleFactor final upperlimit on r is "<<_r95<<endl;
		if(_r95==0) _r95=rmid;

		if(_debug) {start_time=cur_time; cur_time=clock(); cout << "\t\t\tLIMIT_TIMEfor "<<nsize<<" x "<<nexps<<" pseudo exps: " << (cur_time - start_time)/1000000. << " sec\n"; }
		if(ix) delete []ix;
		if(dCLs) delete [] dCLs;
		return _r95;
	}

	double CLsLimit::LimitOnSignalScaleFactor(CountingModel *cms, CLsBase *frequentist, int nexps){
		LimitOnSignalScaleFactor(cms,  1, 1, frequentist, nexps);
	}

	void CLsLimit::DoingStatisticalBandsForCLs(CountingModel *cms, CLsBase *frequentist, int nexps){
		clock_t start_time, cur_time, funcStart_time;
		start_time=clock(); cur_time=clock(); funcStart_time=clock();
		_frequentist = frequentist; _nexps=nexps;

		cms->SetSignalScaleFactor(1.0);
		double cls_mean=0;
		double cp=0; // cumulative p
		_vCLs_Req1.clear(); _vCLs_Req1_CP.clear(); 
		cms->UseAsimovData();
		_frequentist->BuildM2lnQ(cms, _nexps);
		if(_debug){
			start_time=cur_time; cur_time=clock(); cout << "\t TIME_in_BuildM2logQ: "<< _nexps<<" toys, "<< (cur_time - start_time)/1000000./60. << " minutes\n"; fflush(stdout);
		}
		vector<double> vm2logQ_b = _frequentist->Get_m2logQ_b()	;
		vector<double> vm2logQ_sb = _frequentist->Get_m2logQ_sb();		
		nexps=vm2logQ_sb.size();  int nexps_b=vm2logQ_b.size();
		if(nexps!=nexps_b) {
			cout<<"**Error exps for s+b = "<<nexps<<", for b-only="<<nexps_b<<endl;
			return;
		}
		if(nexps<=0){
			cout<<"***Error exps = "<<nexps<<endl;
			return;
		}

		vector<double> qbn, pbn; qbn.clear(); pbn.clear();
		SortAndCumulative(vm2logQ_b, qbn, pbn);	//from small to large	
		vector<double> qsbn, psbn; qsbn.clear(); psbn.clear();
		SortAndCumulative(vm2logQ_sb, qsbn, psbn);		
		int nb = qbn.size();
		int nsb = qsbn.size();
		double *dcls = new double[nb];
		double *dpcls = new double[nb];

		if(_debug){
			start_time=cur_time; cur_time=clock(); cout << "\t TIME_in_SortOut -2lnQ: "<< _nexps<<" toys, "<< (cur_time - start_time)/1000. << " milisec \n"; fflush(stdout);
		}

		int previousStopPoint=0;
		for(int i=0; i<nb; i++){
			double cls=0;//1./(double)_nexps;   // FIXME 
			double pcls=0;
			if(i>0){
				pcls=pbn[i]-pbn[i-1];
			}
			else pcls=pbn[0];

			for(int j=previousStopPoint; j<nsb; j++){
				if(qsbn[j] >= qbn[i]){
					cls= ( (j==0?1:(1-psbn[j-1])) / (i==0?1:(1-pbn[i-1])) );	
					previousStopPoint=j;
					break;
				}
			}
			dcls[i]=cls;
			dpcls[i]=pcls;
		}
		if(_debug){
			start_time=cur_time; cur_time=clock(); cout << "\t TIME_in_CLsBands: "<< _nexps<<" toys, "<< (cur_time - start_time)/1000. << " milisec \n"; fflush(stdout);
		}

		int csize=nb;
		int *i_cls=new int[csize];
		Sort(csize, dcls, i_cls, 0);
		for(int i=0; i<csize; i++){
			_vCLs_Req1.push_back(dcls[i_cls[i]]);
			_vCLs_Req1_CP.push_back( i==0?dpcls[i_cls[0]]:(dpcls[i_cls[i]] + _vCLs_Req1_CP[i-1]) );
			cls_mean+=dcls[i]*dpcls[i];
		}

		if(dcls) delete [] dcls;
		if(dpcls) delete [] dpcls;
		if(i_cls) delete []i_cls;

		if(_debug) {
			int step = _vCLs_Req1_CP.size()/20;
			cout<<"\t In ProjectingCLs, printing out the CLs (r=1) and cummulative p"<<endl;
			cout<<"\t size="<<_vCLs_Req1_CP.size()<<", step="<<step<<endl;
			int i=0;
			for(; i<(int)_vCLs_Req1.size(); i+=(step+1))
				printf("%10d\t CLs %0.6f p %7.6f\n",i, _vCLs_Req1[i], _vCLs_Req1_CP[i]);
			if( i!=_vCLs_Req1_CP.size() || (i==_vCLs_Req1_CP.size() && step!=1) ){
				printf("%10d\t CLs %0.6f p %7.6f\n",_vCLs_Req1_CP.size(), _vCLs_Req1.back(), _vCLs_Req1_CP.back());
			}
		}
		GetBandsByFermiCurveInterpolation(_vCLs_Req1,_vCLs_Req1_CP, _CLsProjected[1], _CLsProjected[3], _CLsProjected[0], _CLsProjected[4]);
		_CLsProjected[5]=cls_mean; _CLsProjected[2]=GetBandByFermiCurveInterpolation(_vCLs_Req1, _vCLs_Req1_CP, 0);
		if(_debug) cout<<" \t ProjectingCLs_m2sigma_p2sigma: -2s= "<<_CLsProjected[0]<<" -1s= "<<_CLsProjected[1]<<" mean= "<<cls_mean<<" 1s= "<<_CLsProjected[3]<<" 2s= "<<_CLsProjected[4]<<endl;
		if(_debug){
			start_time=cur_time; cur_time=clock(); cout << "\t TIME_in_CLsBandsEnd: " << (cur_time - funcStart_time)/1000000./60. << " minutes\n"; fflush(stdout);
		}

	}
	void CLsLimit::DoingStatisticalBandsForLimit(CountingModel *cms, CLsBase *frequentist, int nexps, int npossibleoutcomes){
		_frequentist=frequentist; _nexps=nexps; _npossibleoutcomes=npossibleoutcomes;
		cms->SetSignalScaleFactor(1.);

		// background only ------ ...........
		double r_mean=0;
		double cp=0; // cumulative p
		_vR95_CP.clear(); _vR95.clear(); _differentialR95s.clear();
		vector<double> vr_tmp; vr_tmp.clear();
		vector<double> vp_tmp; vp_tmp.clear();


		if( _debug >= 10 ){ // print out s,b se, be, correlation, err_pdf_type
			cout<<"\t CLsLimit::DoingStatisticalBandsForLimit print out details on signal, background and uncertainties ...."<<endl;
			cms->Print();
		}

		vector< vector<int> > vvPossibleOutcomes; vvPossibleOutcomes.clear();
		vector<int> vEntries; vEntries.clear();  // index of vvPossibleOutcomes are the same as vEntries
		vector<int> vbkg_tmp; 

		int _nOutcomes=_npossibleoutcomes;
		if(_debug) cout<<" _nOutcomes="<<_nOutcomes<<endl;

		for(int n=0; n<_nOutcomes; n++){
			vbkg_tmp.clear(); 
			vbkg_tmp=cms->GetToyData_H0();
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
		}	// _nOutcomes

		// rank the possible outcomes and sort it
		int npsbl_outcomes = vvPossibleOutcomes.size();
		int *Entries_v = new int[npsbl_outcomes];
		for(int i=0; i<npsbl_outcomes; i++){
			Entries_v[i]=vEntries.at(i);
		}
		int *iEntries_v = new int[npsbl_outcomes];
		Sort(npsbl_outcomes, Entries_v, iEntries_v, 1); // sorted it by "down",  Desc
		if(_debug) {
			cout<<"\n\t\t CLsLimit::ProjectingLimits possible outcomes and their entries "<<endl;
			cout<<"\t\t n_possible_outcomes size= "<<npsbl_outcomes<<endl;
			for(int i=0; i<npsbl_outcomes; i++){
				printf("\t\t ");
				for(int j=0; j<vvPossibleOutcomes.at(iEntries_v[i]).size(); j++){
					printf("%4d ", vvPossibleOutcomes.at(iEntries_v[i]).at(j));
				}		
				printf("\t %5d\n", vEntries.at(iEntries_v[i]));
			}
		}

		double nchs=cms->NumOfChannels();
		int tmpcount;
		if(npsbl_outcomes<50) tmpcount = int(npsbl_outcomes/10);
		else {
			tmpcount = int(npsbl_outcomes/100);
		}
		for(int n=0; n<npsbl_outcomes; n++){
			if( tmpcount==0 || (tmpcount!=0 && (n%tmpcount == 0)) ) printf(" ... ... ... process %3.0f\%\n", n/(double)npsbl_outcomes*100);
			double p=vEntries.at(iEntries_v[n])/(double)_nOutcomes;
			if(_debug){
				cout<<"CLsLimit::ProjectingRLimits scanning outcomes: ";
				for(int j=0; j<nchs; j++){
					printf("%4d ", vvPossibleOutcomes.at(iEntries_v[n]).at(j));
				}		
				cout<<"\t p="<<p<<endl;
			}
			cms->SetData( vvPossibleOutcomes.at(iEntries_v[n]) ); 

			double r;

			//--------------calc r95% with (s,b,n) by throwing pseudo experiments and using the fits to get r95% at CLs=5
			r=LimitOnSignalScaleFactor(cms, _frequentist, _nexps);

			if(_debug) {
				vector<double> vvr, vvcls;
				vvr.clear(); vvcls.clear();
				vvr=GetvTestedScaleFactors();	
				vvcls=GetvTestedCLs();
				cout<<" CLsLimit::DoingStatisticalBandsForLimit r95\% for n="<<n<<", printing out the recorded r and cls:"<<endl;
				for(int i=0; i<(int)vvr.size(); i++)
					printf("r %5.2f cls %3.2f\n",vvr[i], vvcls[i]);
			}

			vr_tmp.push_back(r);
			r_mean+=r*p;
			cp+=p;//Poisson(n,b);
			vp_tmp.push_back(p);

			int nentries_for_thisR=(int)(vEntries.at(iEntries_v[n]));
			for(int ntmp=0; ntmp<nentries_for_thisR; ntmp++){
				_differentialR95s.push_back(r);
			}

			if(_debug)cout<<"n="<<n<<" r95="<<r<<" p="<<p<<" cummulative p="<<cp<<endl;
		}

		if(Entries_v) delete [] Entries_v;
		if(iEntries_v) delete [] iEntries_v;

		int nr = vr_tmp.size(); 
		double *r_v = new double[nr];	
		int *ir_v = new int[nr];	
		for(int i=0; i<nr; i++) r_v[i]=vr_tmp.at(i);
		Sort(nr, r_v, ir_v, 0); // from small to large
		for(int i=0; i<nr; i++){

			_vR95.push_back(vr_tmp.at(ir_v[i]));
			if(i==0) _vR95_CP.push_back(vp_tmp.at(ir_v[i]));
			else _vR95_CP.push_back( _vR95_CP.at(i-1) + vp_tmp.at(ir_v[i]) ); 
		}	
		if(r_v) delete [] r_v;
		if(ir_v) delete [] ir_v;

		if(_debug) {
			cout<<" CLsLimit::ProjectingLimits, printing out the r and cummulative p for all possible outcomes:"<<endl;
			for(int i=0; i<(int)_vR95.size(); i++)
				printf("r %5.2f p %5.4f\n",_vR95[i], _vR95_CP[i]);
		}

		GetBandsByFermiCurveInterpolation(_vR95,_vR95_CP, _r95Projected[1], _r95Projected[3], _r95Projected[0], _r95Projected[4]);
		_r95Projected[5]=r_mean; _r95Projected[2]=GetBandByFermiCurveInterpolation(_vR95, _vR95_CP, 0.5);
		if(_debug) {
			cout<<"CLsLimit::ProjectingLimits projecting r at "<<endl; 
			cout<<"    -2 sigma ="<<_r95Projected[0]<<endl;
			cout<<"    -1 sigma ="<<_r95Projected[1]<<endl;
			cout<<"    median   ="<<_r95Projected[2]<<endl;
			cout<<"     1 sigma ="<<_r95Projected[3]<<endl;
			cout<<"     2 sigma ="<<_r95Projected[4]<<endl;
			cout<<"     mean    ="<<_r95Projected[5]<<endl;
		}
	}
	double CLsLimit::Limit_sigma(int nsigma){ // return R at -2sigma, -1sigma, median(50%), 1sigma, 2sigma
		if(nsigma<-2 || nsigma>2) 
		{cout<<"CLsLimit::R_sigma  Error, nsigma should be from -2 to 2, return 0"<<endl;return 0;} 
		else return _r95Projected[nsigma+2];
	}
	double CLsLimit::CLs_sigma(int nsigma){ // return CLs at -2sigma, -1sigma, median(50%), 1sigma, 2sigma
		if(nsigma<-2 || nsigma>2) 
		{cout<<"CLsLimit::CLs_sigma Error, nsigma should be from -2 to 2, return 0"<<endl;return 0;} 
		else return _CLsProjected[nsigma+2];
	}
	void CLsLimit::SetAlpha(double alpha) {_alpha=alpha;} // Confidence Level = 1 - alpha
	double CLsLimit::GetLimit(){return _r95;}
	vector<double> CLsLimit::GetvTestedScaleFactors(){return _vR;} //
	vector<double> CLsLimit::GetvTestedCLs(){return _vCLs;}//	
	double CLsLimit::Limit_mean(){return _r95Projected[5];} //return average value mathmatically.....
	vector<double> CLsLimit::GetDifferentialLimits(){return _differentialR95s;}
	vector<double> CLsLimit::GetvLimits(){return _vR95;} // corresponding to all possible outcomes  ,  cummulative
	vector<double> CLsLimit::GetvLimits_CP(){return _vR95_CP;} // corresponding to all possible outcomes
	double CLsLimit::CLs_mean(){return _CLsProjected[5];} //return average value mathmatically.....
	vector<double> CLsLimit::GetDifferentialCLsReq1(){return _differentialCLs_Req1;}
	vector<double> CLsLimit::GetvCLsReq1(){return _vCLs_Req1;} // corresponding to all possible outcomes
	vector<double> CLsLimit::GetvCLsReq1_CP(){return _vCLs_Req1_CP;} // corresponding to all possible outcomes
	void CLsLimit::SetDebug(int debug){_debug=debug;}
	CLsBase* CLsLimit::GetFrequentist(){return _frequentist;}
	void CLsLimit::SetCLsTolerance(double tolerance){_clstolerance=tolerance;}

};
