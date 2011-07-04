#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TFile.h"
#include "TObject.h"
#include "TTree.h"
#include "TH1D.h"
#include "TMath.h"
#include "TPaveText.h"

#include <iostream>
#include <string>
#include <vector>
#include <iomanip>
#include <math.h>
#include <ctime> // upto second
#include <time.h> // upto micro second

//#include "/data/Projects/cvs/UserCode/mschen/LandS/include/UtilsROOT.h"
//gSystem->Load("/data/Projects/cvs/UserCode/mschen/LandS/test/PlotUtilities_standalone_cc.so");
#include "/data/Projects/cvs/UserCode/mschen/LandS/test/PlotUtilities_standalone.cc"


using namespace std;
using TMath::Sort;
double limit(double &_r95err, vector<double>& _vR, vector<double>& _vCLs, vector<double>& vCLsErr);
bool runCLs(double r, double &cls, double &clserr);
bool GetCLs(TString files, double &cls, double &clserr, bool bGetQdataFromFile=true, double qdata=0);
bool GetM2lnQ(TString file, vector<double> &clsb, vector<double>&clb);
bool GetPValue(vector<double> vclsb, double qdata, double &ret, double &err);
TObject* GetTObject(string filename, string objname);
void StringStrip( std::string & str ); 
void StringSplit( std::vector < std::string > & splitValues, 
		const std::string & str,
		const std::string & delim ); 
double LogLinearInterpolation(double x1, double y1, double x2, double y2, double y);
double LogLinearInterpolationErr(double x1, double y1, double ey1, double x2, double y2, double ey2, double y);

int _debug = 10;

int seed=1234;
double epsilon = 0.001; // _clstolerance
double _alpha = 0.05;
int ntoysForCLsb = 100; 
int ntoysForCLb = 100; // x5 = ntoysForCLsb
int ntoysPerRun = 100;

//TString datacards = "card_shape_H140_*.txt";
//TString datacards = "/data/Projects/cvs/HiggsAnalysis/CombinedLimit/data/benchmarks/gammas/counting-B4-Obs4-Syst50B-gmN.txt";
TString datacards = "card_shape_H140_2mu2e.txt";

void runLHC(){
	gSystem->Exec("rm -f m2lnq_*root");
	gSystem->Exec("rm log1 -f");
	int n =0;
	for(double r = 0.1; r<1; r+=(0.05*r)){
		cout<<r<<" ";
		n++;
	}
	cout<<endl;
	cout<<n<<endl;

	return;

	double err;
	vector<double> vR, vCLs, vCLserr;
	double r95 = limit(err, vR, vCLs, vCLserr);

	cout<<"vR.size="<<vR.size()<<" vCLs.size="<<vCLs.size()<<endl;

	TPaveText *pt = SetTPaveText(0.2, 0.7, 0.3, 0.9);
	//DrawEvolution2D d2d(_vR, _vCLs, "; r ; CLs", (jobname+"_r_vs_cl").Data(), pt);
	DrawEvolution2D d2d(vR, vCLs, "; r ; CLs", "r_vs_cl", pt);
	d2d.draw();
	d2d.save();

	cout<<" 95%CL Upper Limit: "<<r95<<"+/-"<<err<<endl;
}
bool runCLs(double r, double &cls, double &clserr){
	seed++;

	TString command="lands.exe -M Hybrid -tH 100 --testStat LHC --tossToyConvention 1 --UseBestEstimateToCalcQ 0 -v 1 --bNotCalcCLssbb --bSaveM2lnQ --nToysForCLsb 100 --nToysForCLb 100 --fileFittedPars fittedPars.root ";

	//command+=" -d ../hww012jets/input_ee0_mh140.txt ";
	//command+=" -d  card_shape_H140_*.txt ";
	command+=" -d "; command+=datacards;
	command+=" --singlePoint ";
	command+=r; command+=" ";


	TString tmp = command; tmp+=" --bWritePars "; tmp+=" --seed "; tmp+=seed; 
	tmp+=" --fileM2lnQ m2lnq_"; tmp+=seed; tmp+=".root >>log1 ";

	gSystem->Exec(tmp);
	for(int n = 0; n<ntoysForCLsb; n+=ntoysPerRun){
		seed++;
		tmp= command; 
		tmp+=" --bReadPars ";
		tmp+=" --bNotCalcQdata ";
		tmp+=" --seed ";
		tmp+=seed;
		tmp+=" --fileM2lnQ m2lnq_";
		tmp+=seed;
		tmp+=".root >> log1";
		gSystem->Exec(tmp);
	}


	// copy the LimitExtraction Procedure here  ..... ..... ... 
	// still need to chase out why roofit slow  .. ... . . .. . .. . .. ...  . if possible to reduce fitting timing .... .. . .. . .
	GetCLs("m2lnq_*.root", cls, clserr);
	cout<<"CLs: "<<cls<<"+/-"<<clserr<<endl;

	command = "hadd combined_m2lnq_r"; command+=r; command+=".root m2lnq_*root >> log1";
	gSystem->Exec(command);
	command = "rm -f m2lnq_*root";
	gSystem->Exec(command);
	return true;
}
double limit(double &_r95err, vector<double>& _vR, vector<double>& _vCLs, vector<double>& vCLsErr){

	double _r95, cl0, errs0, r1, cl1, errs1;
	//vector<double> _vCLs, _vR, vCLsErr;

	clock_t start_time, cur_time;
	start_time=clock(); cur_time=clock();

	gSystem->Exec("lands.exe -d "+datacards+" -M Bayesian -tB 10 >& logbys");	
	TString bysres = gSystem->GetFromPipe("grep Observed logbys");	
	cout<<bysres<<endl;
	vector<string> tmps;
	StringSplit(tmps, bysres.Data(), " ");
	cout<<tmps.size()<<endl;
	for(int i=0; i<tmps.size(); i++) cout<<tmps[i]<<" "<<endl; 
	cout<<tmps[0]<<endl;
	cout<<"Bayesian Limit = "<<tmps[11]<<endl;

	double r0= TString(tmps[11]).Atof();
	if(r0==0){r0=1;};
	if(_debug){
		start_time=cur_time; cur_time=clock(); cout << "\t TIME_in_BayesianLimitEsitmate: "<< (cur_time - start_time) << " microsec\n";
	}

	if(_debug)cout<<"CLsLimit::LimitOnSignalScaleFactor: start with estimate from Bayesian technique, r95%="<<r0<<endl;

	runCLs(r0, cl0, errs0);

	_vR.push_back(r0);_vCLs.push_back(cl0);
	vCLsErr.push_back(errs0);
	if(_debug)cout<<"Estimated_initial r="<<r0<<"  CLs="<<cl0<< " +/- "<<errs0<<endl;

	//	r1=r0*0.90; //----------usually, CLs-limit is more aggresive than Bayesian's, about 10% smaller.
	if(cl0>_alpha) r1=r0*1.10;	
	else r1=r0*(cl0/_alpha);
	//else r1=r0*0.90;

	runCLs(r1, cl1, errs1);

	_vR.push_back(r1);_vCLs.push_back(cl1);
	vCLsErr.push_back(errs1);
	if(_debug)cout<<"Estimated_r="<<r1<<"  CLs="<<cl1<<" +/- " << errs1<<endl;
	if(fabs(cl1-_alpha)<=epsilon) {
		_r95=r1;
		if(_debug)cout<<"Converge at CLs="<<_alpha<<"+/-"<<epsilon<<" by "<<"2 iterations"<<endl;
		_r95err = LogLinearInterpolationErr(r0, cl0, errs0, r1, cl1, errs1, _alpha);
		return r1;
	}
	if(fabs(cl0-_alpha)<=epsilon) {
		_r95=r0;
		if(_debug)cout<<"Converge at CLs="<<_alpha<<"+/-"<<epsilon<<" by "<<"2 iteration"<<endl;
		_r95err = LogLinearInterpolationErr(r0, cl0, errs0, r1, cl1, errs1, _alpha);
		return r0;
	}

	bool foundit=false;
	double rmid=0;
	int nmaxrepeat=0;	
	while(!foundit && nmaxrepeat<=30 ){ 
		// --- -  --
		//  using linear interpolation to do converge will be quicker
		// ---------
		rmid=LogLinearInterpolation(r0,cl0,r1,cl1,_alpha);

		_r95err = LogLinearInterpolationErr(r0, cl0, errs0, r1, cl1, errs1, _alpha);

		if(rmid<0) rmid=-rmid;  // for CLs limit, constrain r to be > 0,    not for CLsb
		if(_debug >= 10 )cout<<" r0="<<r0<<" cl0="<<cl0<<" r1="<<r1<<" cl1="<<cl1<<" rmid="<<rmid<<endl;
		double clmid, errsmid;
		runCLs(rmid, clmid, errsmid);

		_vR.push_back(rmid);_vCLs.push_back(clmid);
		vCLsErr.push_back(errsmid);
		if(_debug)cout<<"TESTED r="<<rmid<<"  CLs="<<clmid<<" +/- "<<errsmid<<endl;
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
				// there is some improvement on converging in CVS http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/mschen/LandS_diffWaysTossToys/?hideattic=0&sortby=date
				for(int icls= 0; icls<_vR.size(); icls++){
					if(_vCLs[icls]<=_alpha) hasCLsLT05=true;
					if(_vCLs[icls]>=_alpha) hasCLsGT05=true;
				}	
				if(!hasCLsLT05 && _debug) cout<<"\t we don't have CLs < "<<_alpha<<endl;
				if(!hasCLsGT05 && _debug) cout<<"\t we don't have CLs > "<<_alpha<<endl;
				if(!hasCLsLT05) {
					if(rmid>0)rmid *= 1.05; //FIXME this number should more smart 
					if(rmid<0)rmid *= 0.95; //FIXME this number should more smart 

					runCLs(rmid, clmid, errsmid);

					_vR.push_back(rmid);_vCLs.push_back(clmid);
					vCLsErr.push_back(errsmid);
					if(_debug)cout<<"TESTED r="<<rmid<<"  CLs="<<clmid<<" +/- "<<errsmid<<endl;
				}
				if(!hasCLsGT05) {
					if(rmid<0)rmid *= 1.05; //FIXME this number should more smart 
					if(rmid>0)rmid *= 0.95; //FIXME this number should more smart 

					runCLs(rmid, clmid, errsmid);

					_vR.push_back(rmid);_vCLs.push_back(clmid);
					vCLsErr.push_back(errsmid);
					if(_debug)cout<<"TESTED r="<<rmid<<"  CLs="<<clmid<<" +/- "<<errsmid<<endl;
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
		//if(dCLs[ix[i]]==_alpha) {_r95=_vR[ix[i]]; return _r95;}
		if(dCLs[ix[i]]<_alpha) { 
			x1=_vR[ix[i]]; 
			errs0 = vCLsErr[ix[i]];
			y1=dCLs[ix[i]];
		}
		if(dCLs[ix[i]]>=_alpha && x2==0){
			x2=_vR[ix[i]]; 
			errs1 = vCLsErr[ix[i]];
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
	if(_debug) cout<<"\t closest_to_alpha r="<<rmid<<" cls="<<_vCLs.at(tmp_min_index)<<" +/- "<<vCLsErr.at(tmp_min_index)<<endl;

	if(x1==x2 && _debug ) cout<<"\t Warning: r1=r2, we are using the same r to do interpolation, just because we get it twice and one for CLs>alpha, and the other for CLs<alpha"<<endl;

	_r95=LogLinearInterpolation(x1,y1,x2,y2,_alpha);
	_r95err=LogLinearInterpolationErr(x1,y1, errs0, x2,y2, errs1, _alpha);
	if(_debug) cout<<"CLsLimit::LimitOnSignalScaleFactor final upperlimit on r is "<<_r95<< " +/- "<<_r95err<<endl;
	if(_r95==0) _r95=rmid;

	if(_debug) {start_time=cur_time; cur_time=clock(); cout << "\t\t\tLIMIT_TIMEfor "<<nsize<<" x "<<ntoysForCLsb+ntoysForCLb<<" pseudo exps: " << (cur_time - start_time)/1000000. << " sec\n"; }
	if(ix) delete []ix;
	if(dCLs) delete [] dCLs;


	return _r95;

}
bool GetCLs(TString files, double &cls, double &err, bool bGetQdataFromFile, double qdata){
	TString sfiles = gSystem->GetFromPipe("ls "+files);
	cout<<sfiles<<endl;
	vector<string> vsfiles;
	StringSplit(vsfiles, sfiles.Data(), "\n");
	//cout<<vsfiles.size()<<endl;
	vector<double> vclsb, vclb; vclsb.clear(); vclb.clear();
	for(int i=0; i<vsfiles.size(); i++){
		GetM2lnQ(vsfiles[i], vclsb, vclb);
	}
	cout<<" nclsb="<<vclsb.size()<<" nclb="<<vclb.size()<<endl;

	if(bGetQdataFromFile) {
		TH1D *h = (TH1D*)GetTObject(vsfiles[0], "value");
		qdata = h->GetBinContent(2);
	}
	cout<<" qdata = "<<qdata<<endl;

	double clsb, clb, errsb, errb;
	GetPValue(vclsb, qdata, clsb, errsb);
	GetPValue(vclb, qdata, clb, errb);

	if(clb==0){if(_debug)	cout<<"CLsBase::CLs  Warning clb_b==0 !!!!"<<endl; err = 1;  return 1;}
	err = sqrt( errb/clb*errb/clb + errsb/clsb*errsb/clsb) * clsb/clb;
	if(_debug>=10) cout<<"CLsBase::CLs  CLs=CLsb/CLb="<<clsb/clb<<"+/-"<<err<<endl;
	cls = clsb/clb;

	/*
	 * check error derivation
	 TH1D *hsb = new TH1D("hsb", "hsb", 1, 0, 1);
	 TH1D *hb = new TH1D("hb", "hb", 1, 0, 1);
	 hsb->SetBinContent(1, clsb);
	 hb->SetBinContent(1, clb);
	 hsb->SetBinError(1, errsb);
	 hb->SetBinError(1, errb);
	 hsb->Divide(hb);
	 cout<<hsb->GetBinContent(1)<<"+/-"<<hsb->GetBinError(1)<<endl;;
	 */
	return true;
}
bool GetM2lnQ(TString file, vector<double> &vclsb, vector<double>&vclb){
	TTree *tsb = (TTree*)GetTObject(file.Data(), "T1");
	TTree *tb = (TTree*)GetTObject(file.Data(), "T2");
	double clsb, clb;
	TBranch *brCLsb ;
	tsb->SetBranchAddress("brT", &clsb, &brCLsb);
	Long64_t nentries = tsb->GetEntries();
	for(int i=0; i<nentries; i++){
		tsb->GetEntry(i);
		vclsb.push_back(clsb);
	}
	//for(int i=0; i<vclsb.size(); i++) cout<<" "<<vclsb[i]<<" "<<endl;
	TBranch *brCLb ;
	tb->SetBranchAddress("brT", &clb, &brCLb);
	nentries = tb->GetEntries();
	for(int i=0; i<nentries; i++){
		tb->GetEntry(i);
		vclb.push_back(clb);
	}
	//for(int i=0; i<vclb.size(); i++) cout<<" "<<vclb[i]<<" "<<endl;
	return true;		
}

TObject* GetTObject(string filename, string objname){

	// FIXME need to check if filename is exist, and histoname is exist 
	if( gSystem->AccessPathName(filename.c_str())) {cout<<filename<<" couldn't be found"<<endl; exit(0);};
	TFile *f;
	f = (TFile*)gROOT->GetListOfFiles()->FindObject(filename.c_str());
	if(f==NULL) f=new TFile(filename.c_str());
	TObject *h = (TObject*)f->Get(objname.c_str());
	if(!h) {cout<<"object ["<<objname<<"] in file ["<<filename<<"] couldn't be found"<<endl; exit(0);};
	return h;
}

void StringStrip( std::string & str ) {
	//-------------------------------------------------------------------------------
	// Strip spaces and tab at the beginning and the end from a string.
	size_t sPos = 0;
	size_t ePos = str.length();
	while ( str[sPos] == ' ' || str[sPos]=='\t' ) { ++sPos; }
	while ( str[ePos] == ' ' || str[ePos]=='\t' ) { --ePos; }
	str = str.substr( sPos, ePos - sPos );
}

void StringSplit( std::vector < std::string > & splitValues, 
		const std::string & str,
		const std::string & delim ) {
	//-------------------------------------------------------------------------------
	// Split a string by a delimiter and return it's vector of strings.
	std::string str2 = str;

	size_t pos = 0;

	while (( pos = str2.find_first_of( delim )) != std::string::npos ) {
		std::string s = str2.substr(0, pos);
		StringStrip( s );
		if(s!="")splitValues.push_back( s );
		str2 = str2.substr( pos + delim.length());
	}

	StringStrip( str2 );
	if(str2=="")return;
	splitValues.push_back( str2 );
}
bool GetPValue(vector<double> vclsb, double qdata, double &ret, double &err){
	double tmp = qdata;
	int _nexps = vclsb.size();
	for(int i=0; i<_nexps; i++){
		if(vclsb[i]>=tmp)
			ret ++ ;	
	}		
	ret/=_nexps;

	err= sqrt(ret*(1-ret)/_nexps);
	if(ret==0||ret==1) err= 1./_nexps;

	if(_debug>=10){
		cout<<"p ="<<ret<<" +/- "<<err<<" and total exps="<<_nexps<<endl;
		if(ret*_nexps <= 20) cout<<"p*nexps="<<ret*_nexps<<", statistic may not enough"<<endl;
	}
	if(ret == 0){
		if(_debug)	cout<<"p=0, it means number of pseudo experiments is not enough"<<endl;
		if(_debug)	cout<<"              Currently, we put p=1./"<<_nexps<<endl;
		ret = 1./(double)_nexps;
	}
	return true;
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
