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
#include <utility>

#include "CLsLimit.h"
#include "CRandom.h"
#include "Utilities.h"
#include "UtilsROOT.h"
#include "BayesianBase.h"


#include "TMinuit.h"

#include "TMath.h"
#include "TH1D.h"
#include "TFile.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TGraphErrors.h"


#include "RooAbsData.h"
#include "RooWorkspace.h"

using std::cout;
using std::endl;
using std::min;
namespace lands{
	CountingModel *cms_global = 0;
	vector<double> vdata_global;
	TMinuit *myMinuit = 0;
	TF1 * fitToRvsCL_expo = new TF1("fitToRvsCL_expo","[0]*exp([1]*(x-[2]))", 0, 0);
	double _signalScale=0;
	double * _inputNuisances = 0;	
	double * _inputNuisancesSigma = 0;	
	double * _startNuisances = 0;	
	double * _minNuisances = 0; 
	double * _maxNuisances = 0; 
	bool _bPositiveSignalStrength = true;
	vector< vector< vector<float> > > vvv_cachPdfValues;
	vector< vector< double > > vv_cachCountingParts;

	vector<double> _lastParams;
	vector<double> _currParams;

	double _customRMin = 0;
	double _customRMax = 0;

	double _countPdfEvaluation = 0;

	bool _bDumpFinalFitResults = 0;

	double del_oldn = 0; double del_newn =0 ; //DELETEME
	void Chisquare(Int_t &npar, Double_t *gin, Double_t &f,  Double_t *par, Int_t iflag){

		int debug = cms_global->GetDebug();
		// par[0] for the ratio of cross section, common signal strength ....
		if(!cms_global)  {
			cout<<"cms_global not pointed yet "<<endl;
			exit(0);
		}

		// DELETEME
		if( (cms_global->Get_max_uncorrelation()>2) && debug){
			del_oldn = del_newn;
			del_newn = _inputNuisances[2];
			if(del_newn != del_oldn) cout<<"DELETEME input nuisance changed !!!! old="<<del_oldn<<", new="<<del_newn<<endl;
		}

		f = 0; // fabs(cms_global->GetRdm()->Gaus() ) ;
		if(par[0]<0 && _bPositiveSignalStrength) { f = -9e10; return; }

		if(par[0]>9e10) {f=-9e10; return;}
		for(int i=0; i<npar; i++) {if(isnan(par[i]) or isinf(par[i])) f=9e20; return;}
		
		bool bAllChannelsAreFlagged = false;
		npar = cms_global->Get_max_uncorrelation()+1;
		if(_lastParams.size()==0){
			for(int i=0;i<npar;i++)	_lastParams.push_back(par[i]);
			cms_global->FlagAllChannels();
			bAllChannelsAreFlagged = true;
		}else{
			for(int i=0; i<npar; i++) {
				if(par[i]!=_lastParams[i]){
				       	cms_global->FlagChannelsWithParamsUpdated(i);
					_lastParams[i]=par[i];
				}
			}
		}
		
		/*
		for(int i=0;i<npar; i++){
			if(par[i]>_maxNuisances[i] or par[i]<_minNuisances[i]) { f=-99999; return; }
		}
		*/


		const VChannelVSampleVUncertainty &vvv_idcorrl = (cms_global->Get_vvv_idcorrl());
		const VChannelVSampleVUncertainty &vvv_pdftype = (cms_global->Get_vvv_pdftype());
		const VChannelVSampleVUncertaintyVParameter &vvvv_uncpar = cms_global->Get_vvvv_uncpar();

		cms_global->SetSignalScaleFactor(par[0]); // scale the norminal set of signal normalizations 

		const VChannelVSample &vv_sigbks = cms_global->Get_vv_exp_sigbkgs(); 

		const vector<int> &v_pdftype = cms_global->Get_v_pdftype();
		//cout<<" DELETEME :  Chisquare,  v_pdftype 1 = "<<v_pdftype[1]<<endl;; 
		

		const vector<double> &v_GammaN = cms_global->Get_v_GammaN();
		const vector< vector<double> > &v_paramsUnc = cms_global->Get_v_pdfs_floatParamsUnc();

		Double_t chisq = 0;
		int nchs = cms_global->NumOfChannels();

		if(vdata_global.size() != nchs ) {
			cout<<"vdata_global not set correctly"<<endl;
			cout<<"vdata_global.size = "<<vdata_global.size()<<",  model_channels = "<<nchs<<endl;
			exit(0);
		}
		Double_t tc =0, ss=0,  bs = 0;
		int u=0, s=0, c=0;
		double tmp, tmp2,  ran, h, tmprand;
		int ipar;
		int nsigproc = 1;	
		const double *uncpars;
		bool added = false;
		int indexcorrl, pdftype;
		double norminal = 0;
		double normalization = 0;

		if(vv_cachCountingParts.size()==0){
			vv_cachCountingParts.resize(vv_sigbks.size());
			for(int ch=0; ch<vv_sigbks.size(); ch++){
			//cout<<" DELETEME  vv_sigbkg[ch].size = "<<vv_sigbks[ch].size()<<endl;
				vv_cachCountingParts[ch].resize(vv_sigbks[ch].size());
			}
		}

		for(c=0; c<nchs; c++){
			nsigproc = cms_global->GetNSigprocInChannel(c);	
			tc=0; 
			for(s = 0; s<vvv_pdftype[c].size(); s++){
				if(cms_global->Get_vv_statusUpdated()[c][s]){
					bs = vv_sigbks[c][s];	
					if(cms_global->IsUsingSystematicsErrors()){
						if(cms_global->GetMoveUpShapeUncertainties()){
							const vector<int> &shapeuncs = cms_global->GetListOfShapeUncertainties(c, s);
							h=0;
							added = false;
							for(int i = 0; i<shapeuncs.size(); i++){
								indexcorrl = vvv_idcorrl[c][s][shapeuncs[i]];
								pdftype = vvv_pdftype[c][s][shapeuncs[i]];
								ran = par[indexcorrl];
								uncpars  = &(vvvv_uncpar[c][s][shapeuncs[i]][0]);
								switch (pdftype){
									case typeShapeGaussianLinearMorph:
										if(*(uncpars+7) == 1.){
											tmprand = ran; 
											ran*= (*(uncpars+6));
											if(!added) {h+=*(uncpars+2); added=true; norminal = h; normalization = *(uncpars+3); }
											h += max(-ran, 0.) * (*uncpars) + max(ran, 0.) * (*(uncpars+1)) - fabs(ran)*(*(uncpars+2)) ; // uncer params:  down, up, norminal, normlization_of_main_histogram,  uncertainty_down_onNorm, uncertainty_up_onNorm, scale_down_of_gaussian,  siglebin_or_binned
										}

										break;
									case typeShapeGaussianQuadraticMorph:
										if(*(uncpars+7) == 1.){
											tmprand = ran; 
											ran*= (*(uncpars+6));
											if(!added) {h+=*(uncpars+2); added=true; norminal = h;  normalization = *(uncpars+3); }
											if(fabs(ran)<1.){
												h += ran * (ran-1)/2. * (*uncpars) + ran * (ran+1)/2. * (*(uncpars+1)) - ran*ran*(*(uncpars+2)) ; 
											}else  
												h += max(-ran, 0.) * (*uncpars) + max(ran, 0.) * (*(uncpars+1)) - fabs(ran)*(*(uncpars+2)) ; 
										}
										break;
									default:
										break;
								}
							}

							if(added){
								if(h<=0) { h=10e-9;} // cout<<" *h<0* "<<endl; 
								bs = h*normalization*(s<nsigproc?par[0]:1);
								if(isnan(bs)) {cout<<" morphing bs=nan "<<bs<<" c "<<c<<" s="<<s<<" h=" << h<<" normalization="<<normalization<<endl; }
							}
						}
						for(u=0; u<vvv_pdftype[c][s].size(); u++){
							ran = par[(vvv_idcorrl)[c][s][u]];
							uncpars = &(vvvv_uncpar[c][s][u][0]);
							switch (vvv_pdftype[c][s][u]){
								case typeShapeGaussianLinearMorph:
									if(!cms_global->GetMoveUpShapeUncertainties()){
										if(*(uncpars+7) == 1.){
											tmprand = ran; 
											ran*= (*(uncpars+6));
											h = *(uncpars+2) + max(-ran, 0.) * (*uncpars) + max(ran, 0.) * (*(uncpars+1)) - fabs(ran)*(*(uncpars+2)) ; 
											if(h<=0) h=10e-9;
											if( *(uncpars+2)!=0 && bs!=0) {
												bs*=h/(*(uncpars+2));	
											}else if(bs==0) bs = (*(uncpars+3))*h*(s<nsigproc?par[0]:1);
											else { ;}
											ran = tmprand;
										}
									}
									//bs*=pow( (ran>0? *(uncpars+5):*(uncpars+4)) , ran*(*(uncpars+6)) );
									bs*=pow( (ran>0? *(uncpars+5):*(uncpars+4)) , ran>0?ran*(*(uncpars+6)): -ran*(*(uncpars+6)));
									break;
								case typeShapeGaussianQuadraticMorph:
									if(!cms_global->GetMoveUpShapeUncertainties()){
										if(*(uncpars+7) == 1.){
											tmprand = ran; 
											ran*= (*(uncpars+6));
											if(fabs(ran)<1)
												h = *(uncpars+2) + ran * (ran-1)/2. * (*uncpars) + ran * (ran+1)/2. * (*(uncpars+1)) - ran*ran*(*(uncpars+2)) ; 
											else 
												h = *(uncpars+2) + max(-ran, 0.) * (*uncpars) + max(ran, 0.) * (*(uncpars+1)) - fabs(ran)*(*(uncpars+2)) ;
											if(h<=0) h=10e-9;
											if( *(uncpars+2)!=0 && bs!=0) {
												bs*=h/(*(uncpars+2));	
											}else if(bs==0) bs = (*(uncpars+3))*h*(s<nsigproc?par[0]:1);
											else { ;}
											ran = tmprand;
										}
									}
									//cout<<*(uncpars+5)<<" "<<*(uncpars+4)<<endl;
									bs*=pow( (ran>0? *(uncpars+5):*(uncpars+4)) , ran>0?ran*(*(uncpars+6)): -ran*(*(uncpars+6)));
									//bs*=pow( (ran>0? *(uncpars+4):*(uncpars+5)) , ran>0?ran*(*(uncpars+6)): -ran*(*(uncpars+6)));
									break;
								case typeLogNormal:
									bs*=(pow( 1+ (*(uncpars+ (ran>0?1:0))), ran) );
									if(isnan(bs)) {cout<<" typeLogNormal bs=nan "<<bs<<" c "<<c<<" s="<<s<<" uncpars "<<uncpars[0]<<" "<<uncpars[1]<<" ran="<<ran<<endl;}
									break;
								case typeTruncatedGaussian:
									bs*=(1+ (*(uncpars+ (ran>0?1:0)))*ran);
									break;
								case typeGamma:
									ipar = vvv_idcorrl[c][s][u];
									if(*uncpars>0){
										tmp = vv_sigbks[c][s];	
										if(tmp!=0) { bs/=tmp; bs *= ((*uncpars)*ran*(s<nsigproc?par[0]:1)); }
										else bs = (*uncpars)*ran*(s<nsigproc?par[0]:1);
										//cout<<" DELETEME ****  vv_sigbks[c][s]= "<<vv_sigbks[c][s]<<" ran="<<ran<<"  alpha="<<(*uncpars)<<endl;	
									}else{ // if *uncpars, i.e. rho <0, then it's multiplicative gamma
										bs*=(ran/v_GammaN[ipar]);
									}
									break;	
								case typeFlat:
									bs*=(*uncpars + (*(uncpars+1) - *uncpars)*ran);
									break;
								default:
									cout<<"pdf_type = "<<vvv_pdftype[c][s][u]<<" not defined yet"<<endl;
									exit(0);
							}
							if(isnan(bs)) {cout<<" bs=nan "<<bs<<" c "<<cms_global->GetChannelName(c).c_str()<<" s="<<cms_global->GetProcessNames(c)[s].c_str()<<" pdftype="<<vvv_pdftype[c][s][u]<<endl; }
						}
					}
					vv_cachCountingParts[c][s]=bs;
				}else {
					bs=vv_cachCountingParts[c][s];
				}
				tc+=bs;
				if(isnan(tc)) {cout<<" tc=nan "<<"  bs="<<bs<<" c "<<c<<" s="<<s<<endl; }
			}
			if(vdata_global[c]<=0){
				chisq +=( tc - vdata_global[c]);
				//			chisq +=( tc ); //- vdata_global[c]); // to be identical with ATLAS TDR description, for limit only
			}else { 
				if(tc<=0) {f=10e9;cms_global->FlagAllChannels(); return;} // tc < 0, which means non-physical, return f = 10e9
				chisq += (tc-vdata_global[c] - vdata_global[c]*log(tc/vdata_global[c]));
				if(isnan(chisq)) {cout<<" chisq = nan "<<" tc="<<tc<<" data="<<vdata_global[c]<<" in channel: "<<c<<endl;}
			}
			//		else chisq += (tc - vdata_global[c]*log(tc));   // to be identical with ATLAS TDR description, for limit only
		}
		if(isnan(chisq) || chisq>=10e9){ 
			cout<<" DELETEME counting ** ** ** ** chi2="<<chisq<<endl;
			cms_global->FlagAllChannels();
			f=10e9; return; 
		} // checking if it's nan  
		if(cms_global->hasParametricShape()){
			chisq+=	cms_global->EvaluateChi2(par, vvv_cachPdfValues);// use default, norminal sigbkgs for evaluation, not randomized one 		
		}
		if(isnan(chisq) || chisq>=10e9){ 
			cout<<" DELETEME shaping ** ** ** ** chi2="<<chisq<<endl;
			cms_global->FlagAllChannels();
			f=10e9; return; 
		} // checking if it's nan  
		// to be identical with ATLAS TDR description, for limit only
		//http://cdsweb.cern.ch/record/1159618/files/Higgs%20Boson%20%28p1197%29.pdf
		chisq*=2;  //
		if(cms_global->IsUsingSystematicsErrors()){
			// FIXME  when    unc = 0,  then  don't add it 
			double k, tmp;
			for(u=1; u<=cms_global->Get_max_uncorrelation(); u++){
				//normal distribution (mu=0, sigma=1):  1/sqrt(2*PI) * exp(-x^2/2)
				// we are evaluating chisq = -2 * ( ln(Q_H1) - ln(Q_H0) )
				switch (v_pdftype[u]) {
					case typeTruncatedGaussian:
					case typeShapeGaussianQuadraticMorph:
					case typeShapeGaussianLinearMorph:
					case typeLogNormal:
						chisq += pow(par[u]-_inputNuisances[u],2); // make sure if doing lep/tev type and also data fit, then _inputNuisances = norminal set 
						//if(isnan(chisq)) {cout<<"DELETEME chi2=nan _inputNuisances["<<u<<"]="<<_inputNuisances[u]<<endl;}
						break;
					case typeGamma:
						// this is important, one need constraint on the pdf 
						//see wiki,  gamma distribution with theta = 1 :  x^(k-1)*exp(-x)/ (k-1)!
						//k = v_GammaN[u];
						k = _inputNuisances[u];
						if(par[u]<=0) tmp = -par[u];
						else tmp = (k-1)*log(par[u]) - par[u];

						chisq-=tmp*2; // we are evaluating chisq = -2 * ( ln(Q_H1) - ln(Q_H0) )
						break;
					case typeBifurcatedGaussian:
						{
							//Double_t arg = par[u] - v_paramsUnc[u][0]; // x-mean
							Double_t arg = par[u] - _inputNuisances[u]; // x-mean

							Double_t coef(0.0);

							if (arg < 0.0){
								double sigmaL = v_paramsUnc[u][1];
								if (TMath::Abs(sigmaL) > 1e-30) {
									coef = 1./(sigmaL*sigmaL);
								}
							} else {
								double sigmaR = v_paramsUnc[u][2];
								if (TMath::Abs(sigmaR) > 1e-30) {
									coef = 1./(sigmaR*sigmaR);
								}
							}

							chisq+= arg*arg*coef;
							break;
						}
					case typeFlat:
						// do nothing
						break;
					default:
						break;
				}
			}
		}
		// to be identical with ATLAS TDR description, for limit only
		f=chisq;

		if(isnan(f)) f=10e9; // checking if it's nan  

		if( (cms_global->GetPrintParameterFrom() >= 0) && (cms_global->GetPrintParameterTo() >= cms_global->GetPrintParameterFrom()) ){
			if(npar>=cms_global->GetPrintParameterFrom()) 
				printf("PAR ");
			for(int i=cms_global->GetPrintParameterFrom(); i<=cms_global->GetPrintParameterTo(); i++)
				if(i<npar) printf( " %7.4f ", par[i]);
			if(npar>=cms_global->GetPrintParameterFrom()) 
				printf(" f=%10.4f\n", f);
		}

		if(cms_global->GetDebug()>100){
			cout<<"PARS ";
			for(int i=0; i<npar; i++){
				printf(" %.6f ", par[i]);
			}
			cout<<"  "<<f<<endl;
		}

		cms_global->UnFlagAllChannels(bAllChannelsAreFlagged?1:0);

	}

	double MinuitFit(int model, double &r , double &er, double mu /* or ErrorDef for Minos*/, double *pars, bool hasBestFitted, int debug, int *success ){
				//	for(int i=1; i<=cms_global->Get_max_uncorrelation(); i++) {
				//		cout<<"DELETEMEfitstart par "<<i<<" "<<_inputNuisances[i]<<endl;
				//	}

		//cout<<" ********************   MinuitFit **************** "<<endl;


		RooAbsArg::setDirtyInhibit(1);
		if(!(pars && hasBestFitted)){
			_lastParams.clear();	_currParams.clear();
			vvv_cachPdfValues.clear();
			vv_cachCountingParts.clear();
		}

		_signalScale = cms_global->GetSignalScaleFactor();

		int UseMinos = 0;
		if(model == 102) UseMinos = 2; // PL approximation method using Minos ....  with migrad 
		if(model == 101) UseMinos = 1; // PL approximation method using Minos ....  without migrad
		if(model == 1001) UseMinos = 2; // PL approximation method using Minos ....  with migrad 
		if(model == 201 || model==202) UseMinos = 2; // PL approximation method using Minos ....  with migrad   allowing negative mu

		double minuitStep = 0.1;

		int npars = cms_global->Get_max_uncorrelation();
		if(debug>=10)cout<<" MinuitFit npars =  "<<npars<<endl;
		if( !(cms_global->IsUsingSystematicsErrors())) npars=0;
		if(model == 10 && pars){ // fixing  all parameters and get the chi2 
			int tmp;
			double l;
			Chisquare(tmp, 0, l, pars, 0);
			cms_global->SetSignalScaleFactor(_signalScale);
			return l;
		}
		if( (cms_global->IsUsingSystematicsErrors() && npars>0 )  || model ==2 or model==21 or model==101 or model==102 or model==201 or model == 202){

			//FIXME temporarily solution:  when reading a source with all error = 0,  then assign it to be logNormal, error =0,  in UtilsROOT.cc 
			//good solution: redefine npars here, count only sources with definded pdf. 

			//TMinuit *myMinuit = new TMinuit(npars+2);  //initialize TMinuit with a maximum of 5 params
			if(myMinuit) delete myMinuit;
			myMinuit = new TMinuit(npars+2);  //initialize TMinuit with a maximum of 5 params
			myMinuit->SetFCN(Chisquare);

			Double_t arglist[10];
			Int_t ierflg = 0;

			if(debug<100){
				arglist[0]=-1;
				myMinuit -> mnexcm("SET PRINT", arglist, 1, ierflg);
				myMinuit -> mnexcm("SET NOW", arglist, 1, ierflg);

			}
			if(cms_global->Get_minuitPrintLevel() >=0 ){
				arglist[0]=cms_global->Get_minuitPrintLevel();
				myMinuit -> mnexcm("SET PRINT", arglist, 1, ierflg);
			}

			arglist[0] = cms_global->Get_minuitSTRATEGY(); // set to be 0 (was default 1), reduce lots of function calls, ---> speed up by a factor of 6 
			if(debug)cout<<" SET STRATEGY "<<arglist[0]<<endl;
			myMinuit->mnexcm("SET STRATEGY", arglist ,1,ierflg);
			//myMinuit -> mnexcm("SET NOG", arglist, 1, ierflg); // no gradiants required of FCN

			//myMinuit->SetMaxIterations(500); // doesn't seem to do anything in TMinuit

			// Set starting values and step sizes for parameters
			// myMinuit->mnparm(par_index, "par_name", start_value, step_size, lower, higher, ierflg);
			vector<int> v_pdftype = cms_global->Get_v_pdftype();
			vector<double> v_TG_maxUnc = cms_global->Get_v_TruncatedGaussian_maxUnc();
			vector<double> v_GammaN = cms_global->Get_v_GammaN();
			vector< vector<double> > v_paramsUnc = cms_global->Get_v_pdfs_floatParamsUnc();
			vector<string> v_uncname = cms_global->Get_v_uncname();
			vector<bool> v_uncFloatInFit= cms_global->Get_v_uncFloatInFit();
			double maxunc;

			double nuisancesFitRange  =  cms_global->Get_nuisancesRange(); // default is 5 sigma variation allowed

			if(debug) cout<<" nuisancesFitRange = "<<nuisancesFitRange<<" sigma "<<endl;
			for(int i=1; i<=npars; i++){
				if(debug>=10) cout<<" in par "<<i<<endl;
				TString sname=v_uncname[i-1]; 
				if(debug>=10)cout<<" ******* DELETEME in MinuitFit    "<<sname<<":   v_pdftype = "<< v_pdftype[i]<<"  typeGamma = "<<typeGamma<<endl;
				switch (v_pdftype[i]){
					case typeShapeGaussianLinearMorph:
					case typeShapeGaussianQuadraticMorph:
					case typeLogNormal:
						// FIXME need to be smart here ,  when calc significance, if the S > 5 at the end, print out the WARNING message to change the range setting here 
						// or try to get the option of significance from main program 
						// but now [-20, 20] cause some problem in minuit fitting for non-signifcant deviation .... 
						myMinuit->mnparm(i, sname, hasBestFitted?pars[i]:_startNuisances[i], minuitStep, -nuisancesFitRange, nuisancesFitRange, ierflg); // was 5,  causing problem with significance larger than > 7 
						//myMinuit->mnparm(i, sname, hasBestFitted?pars[i]:_startNuisances[i], minuitStep, -10, 10,ierflg); // was 5,  causing problem with significance larger than > 7 
						break;
					case typeTruncatedGaussian :
						maxunc = v_TG_maxUnc[i];	
						if(maxunc>0.2) maxunc = -1./maxunc;
						else maxunc = -5;   // FIXME is hear also need to be extended to -20  ?
						//myMinuit->mnparm(i, sname, hasBestFitted?pars[i]:_inputNuisances[i], minuitStep, maxunc, 20,ierflg); // was 5
						// FIXME need to be smart here ,  when calc significance, if the S > 5 at the end, print out the WARNING message to change the range setting here 
						//myMinuit->mnparm(i, sname, hasBestFitted?pars[i]:_inputNuisances[i], minuitStep, maxunc, 5,ierflg); // was 5
						myMinuit->mnparm(i, sname, hasBestFitted?pars[i]:_startNuisances[i], minuitStep, maxunc, 5,ierflg); // was 5
						break;
					case typeGamma:
						//myMinuit->mnparm(i, sname, v_GammaN[i], 0.5, 0, 100000, ierflg); // FIXME,  could be 100 times the N if N>0,  100 if N==0
						//myMinuit->mnparm(i, sname, hasBestFitted?pars[i]:_inputNuisances[i], minuitStep, 0, (v_GammaN[i]+1)*5, ierflg); // FIXME,  could be 100 times the N if N>0,  100 if N==0
						myMinuit->mnparm(i, sname, hasBestFitted?pars[i]:_startNuisances[i], minuitStep, 0, (v_GammaN[i]+1)*5, ierflg); // FIXME,  could be 100 times the N if N>0,  100 if N==0
						//myMinuit->mnparm(i, sname, hasBestFitted?pars[i]:_inputNuisances[i], 1, 0, v_GammaN[i]+4*sqrt(v_GammaN[i]+1), ierflg); // FIXME,  could be smarter in +/-5sigma stat range
						break;
					case typeBifurcatedGaussian:
						//myMinuit->mnparm(i, sname, hasBestFitted?pars[i]:_inputNuisances[i], minuitStep, v_paramsUnc[i][3], v_paramsUnc[i][4], ierflg  );
						if(debug>=10) cout<<" _startNuisances = "<<_startNuisances[i]<<endl;
						myMinuit->mnparm(i, sname, hasBestFitted?pars[i]:_startNuisances[i], minuitStep, v_paramsUnc[i][3], v_paramsUnc[i][4], ierflg  );
						break;
					case typeFlat:
						//myMinuit->mnparm(i, sname, hasBestFitted?pars[i]:_inputNuisances[i], minuitStep, 0, 1, ierflg  );
						myMinuit->mnparm(i, sname, hasBestFitted?pars[i]:_startNuisances[i], minuitStep, 0, 1, ierflg  );
						break;
					default:
						cout<<"pdftype not yet defined:  "<<v_pdftype[i]<<", npars="<<npars<<", i="<<i<<endl;
						cout<<"**********"<<endl;
						exit(0);
				}
				if(v_uncFloatInFit[i-1]==false)myMinuit->FixParameter(i);
			}

			//if(debug>=10) cout<<"DELETEME in MinuitFit    1"<<endl;

			// through fixing the ratio to determine whether fit for S+B(r=1) or B-only (r=0)   Q_tevatron
			// let the ratio float, then it's Q_atlas
			if(model==1){ // S+B, fix r
				myMinuit->mnparm(0, "ratio", 1, minuitStep, 0, 300, ierflg);
				myMinuit->FixParameter(0);
			}
			else if(model==0 || model==1001){ // B-only, fix r
				myMinuit->mnparm(0, "ratio", 0.0, minuitStep, -1, 300, ierflg);
				myMinuit->FixParameter(0);
			}
			else if(model==2 or model==201 or model == 202){ // S+B,  float r

				double rmin = -100, rmax = 300; 
				if(_customRMax != _customRMin) {rmin=_customRMin; rmax=_customRMax; }
				if(model==2 or model==201)myMinuit->mnparm(0, "ratio", mu, minuitStep, rmin, rmax, ierflg); // andrey's suggestion, alow mu hat < 0   
				if(model==202)myMinuit->mnparm(0, "ratio", r, minuitStep, rmin, rmax, ierflg); // andrey's suggestion, alow mu hat < 0   
				// mu starting point is now configurable via the argument "mu",   when fitting asimov_b, the starting mu should be 0, elsewhere 1
				// mu fitting range maybe need to be configurable via command line. 
				//myMinuit->mnparm(0, "ratio", 1, 0.1, 0, 100, ierflg);  // ATLAS suggestion,   mu hat >=0:   will screw up in case of very downward fluctuation
				_bPositiveSignalStrength = false;
			}
			else if(model==21 || model==101 || model==102){ // S+B,  float r
				double rmin = 0, rmax = 300; 
				if(_customRMax != _customRMin) {rmin=_customRMin; rmax=_customRMax; if (rmin<0) rmin=0;}
				//myMinuit->mnparm(0, "ratio", _startNuisances[0], minuitStep, 0.0, 300, ierflg);  // ATLAS suggestion,   mu hat >=0:   will screw up in case of very downward fluctuation
				myMinuit->mnparm(0, "ratio", 1, minuitStep, rmin, rmax, ierflg);  // ATLAS suggestion,   mu hat >=0:   will screw up in case of very downward fluctuation
				if(model==102) myMinuit->mnparm(0, "ratio", r, minuitStep, rmin, rmax, ierflg); //make starting r configurable
				_bPositiveSignalStrength = true;
			}
			else if(model==3){ // profile mu
				myMinuit->mnparm(0, "ratio", mu, minuitStep, -100, 300, ierflg);
				myMinuit->FixParameter(0);
				_bPositiveSignalStrength = false;
			}
			else if(model==4){ // only floating mu,  not fit for systematics
				myMinuit->mnparm(0, "ratio", mu, minuitStep, -100, 300, ierflg);
				_bPositiveSignalStrength = false;
				for(int i=1; i<=npars; i++) myMinuit->FixParameter(i);
			}
			else if(model==5){ // no profiling at all, i.e. fix all parameters including strength 
				//	myMinuit->mnparm(0, "ratio", mu, 0.1, -100, 300, ierflg);
				//	for(int i=0; i<=npars; i++) myMinuit->FixParameter(i);

				int tmp;
				double l;
				double *par;
				par = new double[npars+1];
				par[0]=mu;
				for(int i=1; i<=npars; i++){
					par[i] = 1.;  // FIXME  why 1 ?    not 0 ?
				}

				Chisquare(tmp, 0, l, par, 0);
				delete []  par;
				cms_global->SetSignalScaleFactor(_signalScale);
				return l;

			}else {
				cout<<"Model not specified correctly:  0-3"<<endl;
				return 0;
			}


			// Setting POIs to be fixed ,  for scanning LL vs. POIs .  need for 2D scan, like mass vs. cross_section 
			if(cms_global->bFixingPOIs()) {
				int npois = 0;
				vector< std::pair<TString, double> > vsPOIsToBeFixed = cms_global->GetPOIsToBeFixed();
				for(int i=0; i<vsPOIsToBeFixed.size(); i++){
					if(vsPOIsToBeFixed[i].first =="signal_strength") {
						myMinuit->mnparm(0, "ratio", vsPOIsToBeFixed[i].second, 0.1, -100, 9e10, ierflg);
						myMinuit->FixParameter(0);
						continue;
					}
					if(cms_global->GetWorkSpaceVaried()->var(vsPOIsToBeFixed[i].first) != NULL){
						// POIs can be only signal_strength or  some parameters in Workspace, e.g. MH 
						for(int j=1; j<=npars; j++){
							TString sname=v_uncname[j-1]; 
							if(vsPOIsToBeFixed[i].first==sname) {
								switch  (v_pdftype[j]){

									case typeBifurcatedGaussian:
										myMinuit->mnparm(j, sname, vsPOIsToBeFixed[i].second, minuitStep, v_paramsUnc[j][3], v_paramsUnc[j][4], ierflg  );
										break;
									case typeFlat:
										if(v_paramsUnc[j][1]==0) break;
										myMinuit->mnparm(j, sname, (vsPOIsToBeFixed[i].second - v_paramsUnc[j][3])/v_paramsUnc[j][1], minuitStep, 0, 1, ierflg  );
										break;
									default:
										break;
								}
								myMinuit->FixParameter(j);
							}
						}
					}
				}
				cms_global->SetbFixingPOIs(false); 
			}			




			//if(debug>=10) cout<<"DELETEME in MinuitFit    2"<<endl;

			arglist[0] = 1;
			if(model==101 || model==102 || model==202 || model==1001)arglist[0] = mu; // ErrorDef for Minos,  just temporaliry using mu ...
			myMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
			// Now ready for minimization step
			arglist[0] = cms_global->Get_maximumFunctionCallsInAFit(); // to be good at minization, need set this number to be 5000 (from experience of hgg+hww+hzz combination)
			if(debug)cout<<" Maximum Function Calls="<<arglist[0]<<endl;
			arglist[1] = cms_global->Get_minuitTolerance();  // tolerance 
			//arglist[1] = 0.009991;
			if(!UseMinos)myMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
			if(UseMinos==2)myMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
			if(debug || ierflg)cout <<" MIGRAD Number of function calls in Minuit: " << myMinuit->fNfcn << endl;
			if(debug || ierflg)cout <<" MIGRAD return errflg = "<<ierflg<<endl;
			if(debug || ierflg)cout <<" MinuitFit("<<model<<")"<<endl;
			if(debug || ierflg) {
				if(model==2 or model==201 )cout<<" starting mu = "<<mu<<endl;
				if(model==102 or model==202)cout<<" starting mu = "<<r<<endl;
			}

			if(myMinuit->fNfcn==0) {
				double *par;
				par = new double[npars+1];
				double pdel, pdelerr; 
				for(int p=0; p<npars+1; p++) {
					myMinuit->GetParameter(p, pdel, pdelerr);
					par[p]=pdel;
				}
				int tmp;
				double l;
				Chisquare(tmp, 0, l, par, 0);
				cms_global->SetSignalScaleFactor(_signalScale);
				return l;
			}
			//	myMinuit->mnexcm("MINI", arglist ,2,ierflg);
			//	myMinuit->mnexcm("IMPROVE", arglist ,2,ierflg);

			if(success) success[0]=ierflg;

			if(UseMinos){
				arglist[0] = cms_global->Get_maximumFunctionCallsInAFit(); // to be good at minization, need set this number to be 5000 (from experience of hgg+hww+hzz combination)
				arglist[1] = 1; // first parameter : signal strength
				int npois = 1;
				for(int p=1; p<cms_global->POIs().size(); p++) { // 0 is always signal strength
					for(int j=1; j<=npars; j++){
						TString sname=v_uncname[j-1]; 
						if(cms_global->POIs()[p].name ==sname) {
							npois += 1;
							arglist[npois]=j+1;
						}
					}
				}
				//arglist[2] = 2;
				//myMinuit->mnexcm("MINOS", arglist , 2, ierflg);
				myMinuit->mnexcm("MINOS", arglist , npois+1, ierflg);
				/*
				   MINOs  [maxcalls] [parno] [parno] ....      parno starts from 1  instead of 0 ... 
				   Causes a Minos error analysis to be performed on the parameters whose numbers are specified. 
				   If none are specified, Minos errors are calculated for all variable parameters. 
				   Minos errors may be expensive to calculate, but are very reliable since they take account of non-linearities
				   in the problem as well as parameter correlations, and are in general asymmetric. The optional argument
				   specifies the (approximate) maximum number of function calls per parameter requested,
				   after which the calculation will be stopped for that parameter.
				 */
				if(debug || ierflg )cout << " MINOS Number of function calls in Minuit: " << myMinuit->fNfcn << endl;
				if(debug || ierflg )cout << " MINOS return errflg = "<<ierflg<<endl;
				if(ierflg){
					cout<<"WARNING: Minos fit fails, try other options"<<endl;
				}
				if(success) success[0]=ierflg;
				if(debug>=10)myMinuit->mnexcm("SHOW COVariance", arglist, 2, ierflg);
				if(debug>=10)myMinuit->mnexcm("SHOW CORrelations", arglist, 2, ierflg);
			}

			myMinuit->GetParameter(0, r, er);
			if(debug)cout<<"DELETEME before calc error,  r="<<r<<"+/-"<<er<<endl;

			// Print results
			Double_t amin,edm,errdef;
			Int_t nvpar,nparx,icstat;
			Double_t errUp, errLow, errParab=0, gcor=0; 
			if(model==101 or model==102 or model == 202){
				myMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
				myMinuit->mnerrs(0, errUp, errLow, errParab, gcor);
				if(errUp==0 and errLow==0) {
					errUp = er;
					errLow = -er;
				}
				if(debug)cout<<"DELETEME mu: errUp="<<errUp<<" errLow="<<errLow<<" errParab="<<errParab<<" gcor="<<gcor<<endl;
				if(cms_global->POIs().size()>0){
					cms_global->setPOI(0, r, errUp, errLow);

					Double_t errUp1, errLow1, errParab1=0, gcor1=0; 
					// POIs can be only signal_strength or  some parameters in Workspace, e.g. MH 
					for(int p=1; p<cms_global->POIs().size(); p++) { // 0 is always signal strength
						for(int j=1; j<=npars; j++){
							TString sname=v_uncname[j-1]; 
							if(cms_global->POIs()[p].name ==sname) {
								double poi, poierr;
								myMinuit->GetParameter(j, poi, poierr);
								myMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
								myMinuit->mnerrs(j, errUp1, errLow1, errParab1, gcor1);
								if(errUp1==0 and errLow1==0) {
									errUp1 = poierr;
									errLow1 = -poierr;
								}
								if(v_pdftype[j]==typeFlat){
									poi = v_paramsUnc[j][3] + v_paramsUnc[j][1] * poi;
									errUp1 = v_paramsUnc[j][1] * errUp1;
									errLow1 = v_paramsUnc[j][1] * errLow1;
								}
								cms_global->setPOI(j, poi, errUp1, errLow1);
								if(debug)cout<<"DELETEME "<<sname<<": errUp="<<errUp1<<" errLow="<<errLow1<<" errParab="<<errParab1<<" gcor="<<gcor1<<endl;
							}
						}
					}
				}

			}

			double l = myMinuit->fAmin;
			myMinuit->GetParameter(0, r, er);
			if(debug)cout<<"DELETEME r="<<r<<"+/-"<<er<<" fMin="<<l<<endl;
			if(errUp==0 and errLow==0) {
				errUp = er;
				errLow = -er;
			}

			if(debug and UseMinos) cout<<" signal_strength :  [ "<<r+errLow<<"  "<<r+errUp<<" ] "<<endl;
			if(model==101 or model==102 or model==202) {
				er=r;
				r+=errUp;   // for upper limit
				er+=errLow; // for lower limit
			}

			if(debug || pars || ierflg){
				_inputNuisancesSigma = cms_global->Get_norminalParsSigma();

				if(debug || ierflg || _bDumpFinalFitResults )printf("  par                 name         fitted_value                      input_value                start_value        dx/s_in,s_out/s_in      (flatParam +/- error)\n");
				for(int i=0; i<=npars; i++){
					double tmp, tmpe;
					myMinuit->GetParameter(i, tmp, tmpe);
					if(debug || ierflg || _bDumpFinalFitResults ) { 
						printf("  par %30s      %.6f +/- %.6f      %.6f +/- %.6f    %.6f     %.2f, %.2f", i>0?v_uncname[i-1].c_str():"signal_strength", tmp, tmpe, _inputNuisances[i], _inputNuisancesSigma[i], _startNuisances[i],  _inputNuisancesSigma[i]==0?0:(tmp-_inputNuisances[i])/_inputNuisancesSigma[i], _inputNuisancesSigma[i]==0?0:tmpe/_inputNuisancesSigma[i]);
						if(i>0 and cms_global->GetWorkSpace()->var(v_uncname[i-1].c_str()) != NULL and v_pdftype[i]==typeFlat) {
							printf("       %.6f +/- %.6f", tmp*v_paramsUnc[i][1]+v_paramsUnc[i][3], tmpe * v_paramsUnc[i][1] );
						}
						printf("\n");
					}
								if(pars && !hasBestFitted)pars[i] = tmp;
					//if(pars)pars[i] = tmp;
				}
			}

				//	for(int i=1; i<=cms_global->Get_max_uncorrelation(); i++) {
				//		cout<<"DELETEMEfitend par "<<i<<" "<<_inputNuisances[i]<<endl;
				//	}
			cms_global->SetSignalScaleFactor(_signalScale);
			RooAbsArg::setDirtyInhibit(0);
			return l;
		}else{
			double par[1];
			int tmp;
			double l;
			if(model==0) par[0]=0;
			else if (model == 1) par[0]=1;
			else if (model == 3) par[0]=mu;
			else {cout<<"model is 2, but going to fix "<<endl; exit(0);}
			Chisquare(tmp, 0, l, par, 0);
			cms_global->SetSignalScaleFactor(_signalScale);
			return l;
		}
		cms_global->SetSignalScaleFactor(_signalScale);
		RooAbsArg::setDirtyInhibit(0);
		return 0.0;
	}	

	bool DoAfit(double mu, vector<double> vdata, vector<RooAbsData*> vrds, double* pars){
		if(cms_global->GetDebug())cout<<"* DoAfit: start with mu= "<<mu<<endl;
		if(!pars) {
			cout<<" pars = 0,  newing "<<endl;
			pars=new double[cms_global->Get_max_uncorrelation()+1];
		}
		vdata_global = vdata;
		cms_global->SetTmpDataForUnbinned(vrds);
		double tmp, tmpr; //double pars[2];
		MinuitFit(3, tmpr, tmpr, mu, pars, 0, cms_global->GetDebug());
		if(cms_global->GetDebug())cout<<"* DoAfit: end   with mu= "<<pars[0]<<endl;
		return true; // should return success or not
	}

	// Class CLsBase
	CLsBase::CLsBase(){
		Q_b = 0;
		Q_sb = 0;
		iq_sb=0;
		iq_b=0;
		_nsig=0; _nbkg=0; _ndat=0;
		_debug=0;
		_rdm=0;
		test_statistics = 1;
		_lognoverb = 0;
	}

	CLsBase::~CLsBase(){
		if(Q_b) delete [] Q_b;
		if(Q_sb) delete [] Q_sb;
		if(iq_b)delete [] iq_b;
		if(iq_sb)delete [] iq_sb;
		if(_lognoverb) delete [] _lognoverb;
		_rdm=0;
	}
	bool CLsBase::BuildM2lnQ(CountingModel *cms, int nexps, int sbANDb_bOnly_sbOnly, bool reUsePreviousToys){
		cms_global = cms;
		_model=cms;
		BuildM2lnQ(nexps, sbANDb_bOnly_sbOnly, reUsePreviousToys);
	}
	bool CLsBase::BuildM2lnQ(int nexps, int sbANDb_bOnly_sbOnly, bool reUsePreviousToys){  // 0 for sbANDb, 1 for bOnly, 2 for sbOnly
		if(test_statistics==1)prepareLogNoverB();
		BuildM2lnQ_data();
		if(sbANDb_bOnly_sbOnly!=2)BuildM2lnQ_b(nexps);
		if(sbANDb_bOnly_sbOnly!=1)BuildM2lnQ_sb(nexps);

		printM2LnQInfo(sbANDb_bOnly_sbOnly);	

		return true;
	}
	void CLsBase::ProcessM2lnQ(){
		Sort(_nexps, Q_b, iq_b, 0); // rank from small to large
		Sort(_nexps, Q_sb, iq_sb, 0);
		if( ( Q_b_data < Q_b[iq_b[0]] || Q_b_data > Q_sb[iq_sb[_nexps-1]] )  && _debug){ 
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
		if(test_statistics==1)
			for(int i=0; i<_nexps; i++){
				tmp.push_back(-2*(Q_sb[i]-_nsig));
			}
		else
			for(int i=0; i<_nexps; i++){
				tmp.push_back(-Q_sb[i]); // back to real -2lnQ
			}
		return tmp;
	} 
	vector<double> CLsBase::Get_m2logQ_b(){
		vector<double> tmp;
		if(test_statistics==1)
			for(int i=0; i<_nexps; i++){
				tmp.push_back(-2*(Q_b[i]-_nsig));
			}
		else
			for(int i=0; i<_nexps; i++){
				tmp.push_back(-Q_b[i]); // back to real -2lnQ
			}
		return tmp;
	} 
	double CLsBase::CLsb(double &err){
		double ret =0;// 1./(double)_nexps;
		double tmp = Q_b_data;
		for(int i=0; i<_nexps; i++){
			if(Q_sb[iq_sb[i]]<=tmp)
				ret = (i+1)/(double)_nexps;	
		}		

		/*
		cout<<" DELETEME "<<endl;
		for(int i=0; i<_nexps; i++){
			cout<<" CLsb "<<i<<" "<<Q_sb[iq_sb[i]]<<endl;
		}
		*/

		err= sqrt(ret*(1-ret)/_nexps);
		if(ret==0||ret==1) err= 1./_nexps;

		if(_debug>=10){
			cout<<"CLsBase::CLsb  CLsb()="<<ret<<" +/- "<<err<<" and total exps="<<_nexps<<endl;
			if(ret*_nexps <= 20) cout<<"CLsBase::CLsb  CLsb*nexps="<<ret*_nexps<<", statistic may not enough"<<endl;
		}
		if(ret == 0){
			if(_debug)	cout<<"CLsBase::CLsb CLsb=0, it means number of pseudo experiments is not enough"<<endl;
			if(_debug)	cout<<"              Currently, we put CLsb=1./"<<_nexps<<endl;
			ret = 1./(double)_nexps;
		}
		return ret;
	}
	double CLsBase::CLs(double &err){
		double errb, errsb;
		double clb=CLb(errb);
		double clsb=CLsb(errsb);
		if(clb==0){if(_debug)	cout<<"CLsBase::CLs  Warning clb_b==0 !!!!"<<endl; err = 1;  return 1;}
		err = sqrt( errb/clb*errb/clb + errsb/clsb*errsb/clsb) * clsb/clb;
		if(_debug>=10) cout<<"CLsBase::CLs  CLs=CLsb/CLb="<<clsb/clb<<"+/-"<<err<<endl;
		return clsb/clb;
	}
	double CLsBase::CLb(double &err){
		return CLb(Q_b_data, err);
	}

	double CLsBase::PValue(double lnq){
		double ret=0;
		bool hasQ_gt_lnq = false;

		if(test_statistics!=6){ // LEP or Tevatron type
			for(int i=0; i<_nexps; i++){ 
				if(Q_b[iq_b[i]] >= lnq)  {
					ret = i/(double)_nexps;	
					hasQ_gt_lnq=true;
					break;
				}
			}		
		}else{ // Standard ProfiledLikelihood ration  
			for(int i=0; i<_nexps; i++){ 
				if(Q_b[iq_b[i]] > lnq)  {
					ret = i/(double)_nexps;	
					break;
				}else{
					hasQ_gt_lnq=true;
				}
			}	
		}
		if(hasQ_gt_lnq==false) {
			ret= test_statistics==6?(1./(double)_nexps):(1-1./(double)_nexps);
			if(_debug or 1) {
				cout<<"********WARNING********"<<endl;
				cout<<" Toys for b-only hypothesis are NOT enough to evaluate the true significance, "<<endl;
				cout<<" Q_b[0]="<<Q_b[iq_b[0]]<<" Q_b["<<_nexps<<"]="<<Q_b[iq_b[_nexps-1]]
					<<", and tested Q="<<lnq<<endl;	
				cout<<" we set PValue to be 1./_nexps = "<<1-ret<<endl;
			}
		}
		return test_statistics==6?ret:(1-ret);
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
		cout<<" This is for evaluating expected mean significance: "<<endl;
		cout<<" Fraction of logQ from S+B hypothesis larger than maximum logQ_b  = "<<fractionGTmaxlnQb<<endl;
		cout<<" I would like this number to be less than 5% "<<endl;
		// if a given lnQ_data larger than maximum logQ_b, it means that toys of b-only is not enough to evaluate significance, 
		// probably we need increase the toy number one or two magnitudes. 

	}
	double CLsBase::CLb(double lnq, double & err){
		double ret =0;// 1./(double)_nexps;
		double tmp = lnq;
		for(int i=0; i<_nexps; i++){ 
			if(Q_b[iq_b[i]]<=tmp)
				ret = (i+1)/(double)_nexps;	
		}		

		/*
		cout<<" DELETEME "<<endl;
		for(int i=0; i<_nexps; i++){
			cout<<" CLb "<<i<<" "<<Q_b[iq_b[i]]<<endl;
		}
		*/

		err= sqrt(ret*(1-ret)/_nexps);
		if(ret==0||ret==1) err= 1./_nexps;

		if(_debug>=10){
			cout<<"CLsBase::CLb  CLb()="<<ret<<"+/-"<<err<<" and total exps="<<_nexps<<endl;
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
		if(ret*_nexps <= 20 && _debug ) cout<<"CLsBase::CLb CLb*nexps="<<ret*_nexps<<", statistic may not enough,  because of data undershoot w.r.t bkg-only hypothesis"<<endl;
		if( (1-ret)*_nexps <= 20 && _debug ) cout<<"CLsBase::CLb  (1-CLb)*nexps="<<(1-ret)*_nexps<<", statistic may not enough"<<endl;
		if(ret == 0 || ret==1){
			if(_debug && ret==0) cout<<"CLsBase::CLb CLb="<<ret<<", it means number of pseudo experiments is not enough"<<endl;
			if(ret==0) {
				if(_debug)cout<<"             data undershoot w.r.t bkg-only hypothesis,   currently, we put CLb=1./"<<_nexps<<endl;
				ret = 1./(double)_nexps;
			}
			if(ret==1) {
				if(_debug)cout<<"            CLb = 1.  --> data overshoot w.r.t bkg-only hypothesis"<<endl;
				//ret = 1 - 1./(double)_nexps;
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

	double CLsBase::Get_m2lnQ_data(){
		if(test_statistics==1)	return -2*(Q_b_data-_nsig);
		else   return -Q_b_data;
	}

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

		// http://en.wikipedia.org/wiki/Normal_distribution
		//sigma erf(n/sqrt(2))   i.e. 1 minus ...    or 1 in ...
		//1 	0.682689492137 	0.317310507863 	3.15148718753
		//2 	0.954499736104 	0.045500263896 	21.9778945081
		//3 	0.997300203937 	0.002699796063 	370.398347380
		//4 	0.999936657516 	0.000063342484 	15,787.192684
		//5 	0.999999426697 	0.000000573303 	1,744,278.331
		//6 	0.999999998027 	0.000000001973 	506,842,375.7


		//	vector<double> vlogQ_sb, vlogQ_sb_prob;
		//	vlogQ_sb.clear(); vlogQ_sb_prob.clear(); 
		//	SortAndCumulative(Q_sb, _nexps, vlogQ_sb, vlogQ_sb_prob, 0);// sort it by  increased order 
		if(ntoys_for_b > 0)  {
			BuildM2lnQ(ntoys_for_b, 1);
		}

		//CheckFractionAtHighEnd(vlogQ_sb, vlogQ_sb_prob);


		double pvalue=PValue(Q_b_data);
		double significance = Significance(pvalue);

		double tmpn = ntoys_for_b*pvalue;  
		double tmpp = ( tmpn - sqrt(tmpn) )/(double)ntoys_for_b;
		double tmpm = ( tmpn + sqrt(tmpn) )/(double)ntoys_for_b;


		if(tmpn<1.8)  tmpp = tmpn/10./(double)ntoys_for_b;

		tmpp = Significance(tmpp);
		tmpm = Significance(tmpm);

		cout<<" p value of data = "<< pvalue << ",  significance = "<< significance << " +"<<tmpp-significance<<" -"<<significance-tmpm<<endl;
		return significance;
	}
	double CLsBase::SignificanceForData(double qdata, vector<double> vq){

		// the input is -2lnQ.  need to translate to 2lnQ
		vector<double> vlnQ_b;
		for(int i=0; i<vq.size(); i++){
			vlnQ_b.push_back(-vq[i]);
		}
		SetLogQ_b(vlnQ_b);
		double pvalue=PValue(qdata);
		double significance = Significance(pvalue);

		int ntoys_for_b = vq.size();
		if(ntoys_for_b<=0) {cout<<" Your input has 0 toy. "<<endl; return 0;}
		double tmpn = ntoys_for_b*pvalue;  
		double tmpp = ( tmpn - sqrt(tmpn) )/(double)ntoys_for_b;
		double tmpm = ( tmpn + sqrt(tmpn) )/(double)ntoys_for_b;


		if(tmpn<1.8)  tmpp = tmpn/10./(double)ntoys_for_b;

		tmpp = Significance(tmpp);
		tmpm = Significance(tmpm);

		cout<<" p value of data = "<< pvalue << ",  significance = "<< significance << " +"<<tmpp-significance<<" -"<<significance-tmpm<<endl;
		return significance;
	}
	/*
	   double CLsBase::SignificanceAnalytically(){

// need a routine to do uncertain number of loops
double prob_data = 0;
double minus_prob_data = 0;
double lnqtmp = 0, p=0;
for(int i=0; i<100; i++){
for(int j=0; j<100; j++)
for(int k=0; k<100; k++){
// reset data numbers 
cms->AddObservedData(0, i);
cms->AddObservedData(1, j);
cms->AddObservedData(2, k);
frequentist.BuildM2lnQ(cms,1);
lnqtmp = frequentist.Get_m2lnQ_data() ;
p = Poisson(i, totalbkg[0]) * Poisson(j, totalbkg[1]) * Poisson(k, totalbkg[2]) ; 
htemp  -> Fill (lnqtmp, p);
if( lnq_data > lnqtmp )
prob_data += p ;
}
}

}
*/

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

void CLsBase::SetTestStatistics(int ts)	{
	if(ts==1 || ts==2 || ts==3 || ts==31 || ts==4 || ts==5 || ts==6)test_statistics = ts; 
	else {
		cout <<"testStat should be 1 for Q_LEP, 2 for Q_TEV or 3 for Q_ATLAS, 31 allowing mu_hat>mu, 4 for only profiling mu, 5 for LHC, 6 for PL "
			<<", your input is not correct: "
			<<ts<<".  ts is set to be default type Q_LEP"
			<<endl; 
	}
}


// Class CLsLimit	
double CLsLimit::LimitOnSignalScaleFactor(CountingModel *cms,
		double minRtoScan, double maxRtoScan,
		CLsBase *frequentist, int nexps, int nstep ){
	cms_global = cms;
	_frequentist=frequentist; _nexps=nexps; 
	_r95err = 0;

	clock_t start_time, cur_time;
	start_time=clock(); cur_time=clock();

	double epsilon = _clstolerance;
	//	int nsigma=5;

	_vR.clear(); _vCLs.clear();

	if(_debug) cout<<"CLsLimit::LimitOnSignalScaleFactor  looking for C.L. 95% Limit on the ratio ----"<<endl;


	double r0,r1;
	double cl0, cl1;

	double errs0, errs1;

	vector<double> vCLsErr; vCLsErr.clear();

	if(cms->GetPhysicsModel()==typeChargedHiggs) {
		if(minRtoScan == maxRtoScan or  minRtoScan<0 or maxRtoScan>1) {
			cout<<"ERROR: for Charged Higgs searches, please provide initial Br range via \" --minRtoScan xxx --maxRtoScan yyy\", where 0 < xxx < yyy < 1 " <<endl; 
			exit(1);
		}
	}

	if(minRtoScan!=maxRtoScan) {

		if(minRtoScan>=maxRtoScan ||(cms->AllowNegativeSignalStrength()==false && minRtoScan <=0) ) {
			cout<<"Error in LimitOnSignalScaleFactor: (minRtoScan="<<minRtoScan<<") >= (maxRtoScan="<<maxRtoScan<<", exit"<<endl;
			cout<<"please make sure maxRtoScan > minRtoScan "<<endl;
			if(cms->AllowNegativeSignalStrength()==false && minRtoScan<=0) 
				cout<<"please make sure minRtoScan > 0 or SetAllowNegativeSignalStrength(true) for your model"<<endl;
			exit(0);
		}
		if(nstep<1)  {cout<<" steps in autoscan should not less than 1, exit"<<endl; exit(0);}
		if(_debug) cout<<"\t First auto scaning R from  "<<minRtoScan<<" to "<<maxRtoScan<<" in "<<nstep<<" steps"<<endl;

		for(double rmid=minRtoScan; rmid<=maxRtoScan; rmid+=(maxRtoScan-minRtoScan)/(double)nstep){
			cms->SetSignalScaleFactor(rmid);
			if(cms->GetTossToyConvention()==1) {
				DoAfit(rmid, cms->Get_v_data_real(), cms->Get_v_pdfs_roodataset_real(), cms->Get_fittedParsInData_sb());
				cms->Set_fittedParsInData_sb(cms->Get_fittedParsInData_sb());
			}
			frequentist->BuildM2lnQ(cms, nexps); 
			if(_rule == 1)
				cl0=frequentist->CLs(errs0); //_nsigma(nsigma);
			else 
				cl0=frequentist->CLsb(errs0); //_nsigma(nsigma);
			_vR.push_back(rmid);_vCLs.push_back(cl0);
			vCLsErr.push_back(errs0);

			if(_debug)cout<<"TESTED r="<<rmid<<"  CLs="<<cl0<<" +/- "<<errs0<<endl;
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
				errs0 = vCLsErr[ix[i]];
			}
			if(dCLs[ix[i]]>_alpha && x2==0){
				x2=_vR[ix[i]]; 
				y2=dCLs[ix[i]];
				errs1 = vCLsErr[ix[i]];
			}
		}
		if(!x1 || !x2 || !y1 || !y2){
			cout<<"Warning: Your initial estimated R range is not suitable, all CLs values I got are in one side of "<<_alpha<<endl;
			x1 = _vR[0]; x2=_vR.back();
			y1 = _vCLs[0]; y2=_vCLs.back();
			errs0 = vCLsErr[0]; errs1 = vCLsErr.back();
		}
		if(ix) delete []ix;
		if(dCLs) delete [] dCLs;
		r0=x1; r1=x2; cl0=y1; cl1=y2;
	} else {
		cms->SetSignalScaleFactor(1.);
		if(cms->GetTossToyConvention()==1) {
			DoAfit(1, cms->Get_v_data_real(), cms->Get_v_pdfs_roodataset_real(), cms->Get_fittedParsInData_sb());
			cms->Set_fittedParsInData_sb(cms->Get_fittedParsInData_sb());
		}
		BayesianBase bys(cms, 0.05, 1.e-2);
		bys.SetNumToys(100);
		bys.SetDebug(_debug);
		r0= bys.Limit();
		if(r0==0){r0=1;};
		if(_debug){
			start_time=cur_time; cur_time=clock(); cout << "\t TIME_in_BayesianLimitEsitmate: "<< (cur_time - start_time) << " microsec\n";
		}

		if(_debug)cout<<"CLsLimit::LimitOnSignalScaleFactor: start with estimate from Bayesian technique, r95%="<<r0<<endl;

		cms->SetSignalScaleFactor(r0);
		if(cms->GetTossToyConvention()==1) {
			cout<<"Before refit: r=1"<<endl;
			for(int i=0; i<cms->Get_max_uncorrelation(); i++) cout<<"par "<<i<<"     "<<cms->Get_fittedParsInData_sb()[i]<<endl;
			DoAfit(r0, cms->Get_v_data_real(), cms->Get_v_pdfs_roodataset_real(), cms->Get_fittedParsInData_sb());
			cms->Set_fittedParsInData_sb(cms->Get_fittedParsInData_sb());
			cout<<"After refit: r="<<r0<<endl;
			for(int i=0; i<cms->Get_max_uncorrelation(); i++) cout<<"par "<<i<<"     "<<cms->Get_fittedParsInData_sb()[i]<<endl;
		}
		frequentist->BuildM2lnQ(cms, nexps); 
		if(_rule == 1)
			cl0=frequentist->CLs(errs0); //_nsigma(nsigma);
		else 
			cl0=frequentist->CLsb(errs0); //_nsigma(nsigma);
		_vR.push_back(r0);_vCLs.push_back(cl0);
		vCLsErr.push_back(errs0);
		if(_debug)cout<<"Estimated_initial r="<<r0<<"  CLs="<<cl0<< " +/- "<<errs0<<endl;

		//	r1=r0*0.90; //----------usually, CLs-limit is more aggresive than Bayesian's, about 10% smaller.
		if(cl0>_alpha) r1=r0*1.10;	
		else r1=r0*(cl0/_alpha);
		//else r1=r0*0.90;

		cms->SetSignalScaleFactor(r1);
		if(cms->GetTossToyConvention()==1) {
			DoAfit(r1, cms->Get_v_data_real(), cms->Get_v_pdfs_roodataset_real(), cms->Get_fittedParsInData_sb());
			cms->Set_fittedParsInData_sb(cms->Get_fittedParsInData_sb());
		}
		frequentist->BuildM2lnQ(cms, nexps); 
		if(_rule == 1)
			cl1=frequentist->CLs(errs1); //_nsigma(nsigma);
		else 
			cl1=frequentist->CLsb(errs1); //_nsigma(nsigma);
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
	}

	bool foundit=false;
	double rmid=0;
	int nmaxrepeat=0;	
	while(!foundit && nmaxrepeat<=30 ){ 
		// --- -  --
		//  using linear interpolation to do converge will be quicker
		// ---------
		rmid=LogLinearInterpolation(r0,cl0,r1,cl1,_alpha);

		if(cms->GetPhysicsModel()==typeChargedHiggs){
			// ****************** for charged higgs
			if(rmid<0) rmid=0;
			if(rmid>1) rmid=1;
			// ****************** end for charged higgs
		}

		_r95err = LogLinearInterpolationErr(r0, cl0, errs0, r1, cl1, errs1, _alpha);

		if(_rule ==1 && rmid<0) rmid=-rmid;  // for CLs limit, constrain r to be > 0,    not for CLsb
		if(_debug >= 10 )cout<<" r0="<<r0<<" cl0="<<cl0<<" r1="<<r1<<" cl1="<<cl1<<" rmid="<<rmid<<endl;
		cms->SetSignalScaleFactor(rmid);
		rmid = cms->GetSignalScaleFactor(); // if not allow negative r, then the scale factor will not be modified in SetSignalScaleFactor. 
		if(cms->GetTossToyConvention()==1) {
			DoAfit(rmid, cms->Get_v_data_real(), cms->Get_v_pdfs_roodataset_real(), cms->Get_fittedParsInData_sb());
			cms->Set_fittedParsInData_sb(cms->Get_fittedParsInData_sb());
		}
		frequentist->BuildM2lnQ(cms, nexps); 
		double clmid, errsmid;
		if(_rule == 1)
			clmid=frequentist->CLs(errsmid); //_nsigma(nsigma);
		else 
			clmid=frequentist->CLsb(errsmid); //_nsigma(nsigma);
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

					if(cms->GetPhysicsModel()==typeChargedHiggs){
						// ****************** for charged higgs
						if(rmid<0) rmid=0;
						if(rmid>1) rmid=1;
						// ****************** end for charged higgs
					}

					cms->SetSignalScaleFactor(rmid);
					if(cms->GetTossToyConvention()==1) {
						DoAfit(rmid, cms->Get_v_data_real(), cms->Get_v_pdfs_roodataset_real(), cms->Get_fittedParsInData_sb());
						cms->Set_fittedParsInData_sb(cms->Get_fittedParsInData_sb());
					}
					frequentist->BuildM2lnQ(cms, nexps); 
					double clmid;
					if(_rule == 1)
						clmid=frequentist->CLs(errsmid); //_nsigma(nsigma);
					else 
						clmid=frequentist->CLsb(errsmid); //_nsigma(nsigma);
					_vR.push_back(rmid);_vCLs.push_back(clmid);
					vCLsErr.push_back(errsmid);
					if(_debug)cout<<"TESTED r="<<rmid<<"  CLs="<<clmid<<" +/- "<<errsmid<<endl;
				}
				if(!hasCLsGT05) {
					if(rmid<0)rmid *= 1.05; //FIXME this number should more smart 
					if(rmid>0)rmid *= 0.95; //FIXME this number should more smart 

					if(cms->GetPhysicsModel()==typeChargedHiggs){
						// ****************** for charged higgs
						if(rmid<0) rmid=0;
						if(rmid>1) rmid=1;
						// ****************** end for charged higgs
					}

					cms->SetSignalScaleFactor(rmid);
					if(cms->GetTossToyConvention()==1){
						DoAfit(rmid, cms->Get_v_data_real(), cms->Get_v_pdfs_roodataset_real(), cms->Get_fittedParsInData_sb());
						cms->Set_fittedParsInData_sb(cms->Get_fittedParsInData_sb());
					}
					frequentist->BuildM2lnQ(cms, nexps); 
					double clmid;
					if(_rule == 1)
						clmid=frequentist->CLs(errsmid); //_nsigma(nsigma);
					else 
						clmid=frequentist->CLsb(errsmid); //_nsigma(nsigma);
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

	if(_debug) {start_time=cur_time; cur_time=clock(); cout << "\t\t\tLIMIT_TIMEfor "<<nsize<<" x "<<nexps<<" pseudo exps: " << (cur_time - start_time)/1000000. << " sec\n"; }
	if(ix) delete []ix;
	if(dCLs) delete [] dCLs;
	return _r95;
}

double CLsLimit::LimitOnSignalScaleFactor(CountingModel *cms, CLsBase *frequentist, int nexps){
	LimitOnSignalScaleFactor(cms,  1, 1, frequentist, nexps);
}

void CLsLimit::DoingStatisticalBandsForCLs(vector<double> vm2logQ_sb, vector<double> vm2logQ_b){
	clock_t start_time, cur_time, funcStart_time;
	start_time=clock(); cur_time=clock(); funcStart_time=clock();

	double cls_mean=0;
	double cp=0; // cumulative p
	_vCLs_Req1.clear(); _vCLs_Req1_CP.clear(); 
	int nexps=vm2logQ_sb.size();  int nexps_b=vm2logQ_b.size();

	if(nexps<=10 || nexps_b<=10){
		cout<<"***Error: too few toys, vsb.size= "<<nexps<<",  vb.size="<<nexps_b<<endl;
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
		start_time=cur_time; cur_time=clock(); cout << "\t TIME_in_SortOut -2lnQ: b"<< nexps_b <<", sb"<<nexps<<" toys, "<< (cur_time - start_time)/1000. << " milisec \n"; fflush(stdout);
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
		start_time=cur_time; cur_time=clock(); cout << "\t TIME_in_CLsBands -2lnQ: b"<< nexps_b <<", sb"<<nexps<<" toys, "<< (cur_time - start_time)/1000. << " milisec \n"; fflush(stdout);
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
		int step = _vCLs_Req1_CP.size()/20+1;
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
	_CLsProjected[5]=cls_mean; _CLsProjected[2]=GetBandByFermiCurveInterpolation(_vCLs_Req1, _vCLs_Req1_CP, 0.5);
	cout<<" \t ProjectingCLs_m2sigma_p2sigma: -2s= "<<_CLsProjected[0]<<" -1s= "<<_CLsProjected[1]<<" mean= "<<cls_mean<<" 1s= "<<_CLsProjected[3]<<" 2s= "<<_CLsProjected[4]<<endl;
}
// this is deprecated by  the "LimitBands" class
void CLsLimit::DoingStatisticalBandsForLimit(CountingModel *cms, CLsBase *frequentist, int nexps, int npossibleoutcomes){
	_frequentist=frequentist; _nexps=nexps; _npossibleoutcomes=npossibleoutcomes;

	cms_global = cms;

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

	vector< VIChannel > vvPossibleOutcomes; vvPossibleOutcomes.clear();
	vector<int> vEntries; vEntries.clear();  // index of vvPossibleOutcomes are the same as vEntries
	VIChannel vbkg_tmp; 

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
const vector<double>& CLsLimit::GetvTestedScaleFactors(){return _vR;} //
const vector<double>& CLsLimit::GetvTestedCLs(){return _vCLs;}//	
double CLsLimit::Limit_mean(){return _r95Projected[5];} //return average value mathmatically.....
const vector<double>& CLsLimit::GetDifferentialLimits(){return _differentialR95s;}
const vector<double>& CLsLimit::GetvLimits(){return _vR95;} // corresponding to all possible outcomes  ,  cummulative
const vector<double>& CLsLimit::GetvLimits_CP(){return _vR95_CP;} // corresponding to all possible outcomes
double CLsLimit::CLs_mean(){return _CLsProjected[5];} //return average value mathmatically.....
const vector<double>& CLsLimit::GetDifferentialCLsReq1(){return _differentialCLs_Req1;}
const vector<double>& CLsLimit::GetvCLsReq1(){return _vCLs_Req1;} // corresponding to all possible outcomes
const vector<double>& CLsLimit::GetvCLsReq1_CP(){return _vCLs_Req1_CP;} // corresponding to all possible outcomes
void CLsLimit::SetDebug(int debug){_debug=debug;}
CLsBase* CLsLimit::GetFrequentist(){return _frequentist;}
void CLsLimit::SetCLsTolerance(double tolerance){_clstolerance=tolerance;}
void CLsLimit::SetRule(int rule){
	_rule = rule; 
	if(rule!=1 && rule!=2 && rule!=3) { 
		cout<<"frequentist rule should be 1 for CLs, or 2 for CLsb, or 3 for FC.  Your input "<<rule<<" is not defined!"<<endl;
		exit(0);
	}
}

// Class CLsLimit	
double CLsLimit::FeldmanCousins(CountingModel *cms,
		double minRtoScan, double maxRtoScan,
		CLsBase *frequentist, int nexps, int nstep ){

	// test statistics muct be       L(n, mu)/L(n, mu_hat)  or   L(n, mu, theta_hat)/L(n, mu_hat, theta_hat)

	cms_global = cms;
	_frequentist=frequentist; _nexps=nexps; 

	double fAdditionalNToysFactor = 2.;
	bool bAdaptiveSampling = true;

	clock_t start_time, cur_time;
	start_time=clock(); cur_time=clock();


	if(_debug) cout<<"CLsLimit::FeldmanCousins  looking for C.L. 95% Limit on the ratio ----"<<endl;


	if(minRtoScan>=maxRtoScan ||(cms->AllowNegativeSignalStrength()==false && minRtoScan <=0) ) {
		cout<<"Error in FeldmanCousins: (minRtoScan="<<minRtoScan<<") >= (maxRtoScan="<<maxRtoScan<<", exit"<<endl;
		cout<<"please make sure maxRtoScan > minRtoScan "<<endl;
		if(cms->AllowNegativeSignalStrength()==false && minRtoScan<=0) 
			cout<<"please make sure minRtoScan > 0 or SetAllowNegativeSignalStrength(true) for your model"<<endl;
		exit(0);
	}
	if(nstep<1)  {cout<<" steps in autoscan should not less than 1, exit"<<endl; exit(0);}
	if(_debug) cout<<"\t First auto scaning R from  "<<minRtoScan<<" to "<<maxRtoScan<<" in "<<nstep<<" steps"<<endl;

	_FCconstruction.clear();
	vector<double> qs; 
	for(double rmid=minRtoScan; rmid<=maxRtoScan; rmid+=(maxRtoScan-minRtoScan)/(double)nstep){
		cms->SetSignalScaleFactor(rmid);
		double thisTestStatistic = 0;
		double sigma;
		double upperEdgeOfAcceptance, upperEdgeMinusSigma, upperEdgePlusSigma;
		double lowerEdgeOfAcceptance, lowerEdgeMinusSigma, lowerEdgePlusSigma;
		if(bAdaptiveSampling){
			// This adaptive sampling algorithm is imported from RooStats http://root.cern.ch/root/html/src/RooStats__NeymanConstruction.cxx.html
			// the adaptive sampling algorithm wants at least one toy event to be outside
			// of the requested pvalue including the sampling variaton.  That leads to an equation
			// N-1 = (1-alpha)N + Z sqrt(N - (1-alpha)N) // for upper limit and
			// 1   = alpha N - Z sqrt(alpha N)  // for lower limit 
			// 
			// solving for N gives:
			// N = 1/alpha * [3/2 + sqrt(5)] for Z = 1 (which is used currently)
			// thus, a good guess for the first iteration of events is N=3.73/alpha~4/alpha
			// should replace alpha here by smaller tail probability: eg. alpha*Min(leftsideFrac, 1.-leftsideFrac)
			// totalMC will be incremented by 2 before first call, so initiated it at half the value
			int totalMC = int(2./_alpha);
			// user control
			double tmc = double(totalMC)*fAdditionalNToysFactor;
			totalMC = (int) tmc; 
			int additionalMC=0;
			bool bUsePreviousToys = false;

			do{
				// this will be executed first, then while conditioned checked
				// as an exit condition for the loop.

				// the next line is where most of the time will be spent 
				// generating the sampling dist of the test statistic.
				additionalMC = 2*totalMC; //grow by a factor of 2

				frequentist->BuildM2lnQ(cms, additionalMC+ (bUsePreviousToys?totalMC:0), 2, bUsePreviousToys); // 2 for s+b hypothesis only ...
				thisTestStatistic=frequentist->Get_m2lnQ_data();

				totalMC = frequentist->GetNexps();

				if(!bUsePreviousToys) bUsePreviousToys=true;

				qs = frequentist->Get_m2logQ_sb();

				sigma = 1;
				upperEdgeOfAcceptance = InverseCDF(qs, _alpha, 1. - _alpha, sigma, upperEdgePlusSigma);
				sigma = -1;
				InverseCDF(qs, _alpha, 1. - _alpha , sigma, upperEdgeMinusSigma);

				sigma = 1;
				lowerEdgeOfAcceptance = InverseCDF(qs, _alpha, 0, sigma, lowerEdgePlusSigma);
				sigma = -1;
				InverseCDF(qs, _alpha, 0, sigma, lowerEdgeMinusSigma);

				if(_debug) cout << " NeymanConstruction: "
					<< "total MC = " << totalMC <<endl; 
				if(_debug>=10)	cout<< "   this test stat = " << thisTestStatistic << endl
					<< " upper edge -1sigma = " << upperEdgeMinusSigma
						<< ", upperEdge = "<<upperEdgeOfAcceptance
						<< ", upper edge +1sigma = " << upperEdgePlusSigma << endl
						<< " lower edge -1sigma = " << lowerEdgeMinusSigma
						<< ", lowerEdge = "<<lowerEdgeOfAcceptance
						<< ", lower edge +1sigma = " << lowerEdgePlusSigma << endl;
			}while(
					( 
					 (thisTestStatistic <= upperEdgeOfAcceptance &&
					  thisTestStatistic > upperEdgeMinusSigma)
					 || (thisTestStatistic >= upperEdgeOfAcceptance &&
						 thisTestStatistic < upperEdgePlusSigma)
					 || (thisTestStatistic <= lowerEdgeOfAcceptance &&
						 thisTestStatistic > lowerEdgeMinusSigma)
					 || (thisTestStatistic >= lowerEdgeOfAcceptance &&
						 thisTestStatistic < lowerEdgePlusSigma) 
					) && (totalMC < 100./_alpha)
			      );

		}else{
			frequentist->BuildM2lnQ(cms, nexps, 2); // 2 for s+b hypothesis only ...
			qs = frequentist->Get_m2logQ_sb();
			// using default comparison (operator <):
			//sort(qs.begin(), qs.end());

			sigma = 1;
			upperEdgeOfAcceptance = InverseCDF(qs, _alpha, 1. - _alpha, sigma, upperEdgePlusSigma);
			sigma = -1;
			InverseCDF(qs, _alpha, 1. - _alpha , sigma, upperEdgeMinusSigma);

			sigma = 1;
			lowerEdgeOfAcceptance = InverseCDF(qs, _alpha, 0, sigma, lowerEdgePlusSigma);
			sigma = -1;
			InverseCDF(qs, _alpha, 0, sigma, lowerEdgeMinusSigma);

		}

		double q_up;
		/* 
		   q_up = qs.back();
		// 1.   quantile with step function 
		//q_up = qs[int(0.95*nexps)];

		// 2. according to FC paper,  over coverage unavoidable due to discreteness 
		vector<double> qn, pn;
		SortAndCumulative(qs, qn, pn);
		for(int i=0; i<pn.size(); i++){
		if(pn[i]>=0.95) 
		//if(pn[i]>0.95) 
		{
		q_up = qn[(i==pn.size()-1)?i:i+1];
		//q_up = qn[i];
		break;
		}
		}


		cout<<endl<<" p: " ;
		for(int i=0; i<pn.size(); i++){
		cout<<pn[i]<<" " ;
		}
		cout<<endl;

		cout<<"p: ";
		for(int i=0; i<20; i++){
		cout<<TMath::Poisson(i, rmid+ vdata_global[0])<<" ";
		}
		cout<<endl;
		*/	

		q_up = upperEdgeOfAcceptance;
		if(_debug){
			cout<<"TESTED r="<<rmid<<"  -2lnQ_up= "<<q_up<<" -2lnQ_data="<< frequentist->Get_m2lnQ_data();
			if(q_up==frequentist->Get_m2lnQ_data())cout<<" =  "<<endl;
			if(q_up>frequentist->Get_m2lnQ_data())cout<<"  >  "<<endl;
			if(q_up<frequentist->Get_m2lnQ_data())cout<<" <  "<<endl;
			//		cout<<" ------  rightmost "<<qn.back()<<endl;
			//		for(int i=0; i<qn.size(); i++){
			//			cout<<" "<<qn[i];
			//		}
			//		cout<<endl;
		}

		qs.push_back(rmid);
		qs.push_back(q_up);
		qs.push_back(frequentist->Get_m2lnQ_data());
		_FCconstruction.push_back(qs);

	}
	// extract the 95% CL  upper limit
	for(int i=0; i<_FCconstruction.size(); i++){
		int j = _FCconstruction.size() - i -1;	
		int n = _FCconstruction[j].size();
		//if(_FCconstruction[j][n-2] > _FCconstruction[j][n-1]) 
		if(_FCconstruction[j][n-2] >= _FCconstruction[j][n-1]) 
		{
			_r95 = _FCconstruction[j][n-3];
			break;
		}
	}

	return _r95;
}

bool CLsBase::BuildM2lnQ_b(int nexps, bool reUsePreviousToys, bool bWriteToys){  // 0 for sbANDb, 1 for bOnly, 2 for sbOnly, 3 for data only 
	RooWorkspace * wtmp = new RooWorkspace("w");
	_inputNuisances = _model->Get_norminalPars();	
	_startNuisances= _model->Get_norminalPars();	
	// effort for adaptive sampling
	int oldNexps = _nexps;// need to be changed   to two parts :   sb toys,  b toys 
	vector<double> tmpQb;
	if(reUsePreviousToys){
		// you have to make sure in the same model with same signal scale factor ...
		if(!Q_b) reUsePreviousToys = false;  // the previous toys are either not exist or deleted
		if(nexps<=oldNexps) reUsePreviousToys = false; // if the new total nexps required is less than previous number ... 
		if(reUsePreviousToys){
			tmpQb.clear();
			for(int i=0; i<oldNexps; i++){
				tmpQb.push_back(Q_b[i]);
			}
		}
	}


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
	_nchannels = _model->NumOfChannels(); // only channels with counting exps,   no parametric shape channel here

	if(Q_b) delete [] Q_b;
	if(iq_b) delete [] iq_b;
	Q_b=new double[_nexps];
	iq_b = new int[_nexps];	

	checkFittedParsInData();

	int nsbi, nbi;
	int tenth = _nexps/10;
	if(tenth<1) tenth=1;
	int ntemp = _nexps*_nchannels;
	if(_debug >=10 ) cout<<"ntemp="<<ntemp<<endl;
	if(test_statistics!=1 && test_statistics!=4){
		ntemp *= _model->Get_max_uncorrelation(); // if using Q_tev or Q_atlas, then multiply by the number of nuisance parameters
		ntemp *= 100;
	}
	if( ntemp>=10000000 ) {
		cout<<"\t gonna generate "<<ntemp*2<<" poisson numbers "<<endl;
	}

	clock_t toytime_start = clock(); int ntoysFor10sec = _nexps;
	for(int i=0; i<_nexps; i++){
		if(i==1) {
			clock_t toytime_stop = clock(); int timeForOneToy = toytime_stop - toytime_start; 
			if(_debug) cout<<" time per toy = "<<timeForOneToy<<" microsec"<<endl;
			if(timeForOneToy>0) ntoysFor10sec = 10*1000000/timeForOneToy+1; 
		}
		if( ntemp>=10000000 or _debug) {
			if( (i+1)%tenth == 0 ){
				printf("... Building -2lnQ,  %4.1f \%\n", i/(double)_nexps*100);
				fflush(stdout);
			}
		}
		Q_b[i]=0;	
		if(reUsePreviousToys && i<oldNexps){
			Q_b[i] = tmpQb[i];
			continue;
		}
		if(_debug){
			if( (i+1)%ntoysFor10sec== 0 ){
				clock_t toytime_stop = clock();
				cout<< " from_1st toy to "<<i+1<<" toy takes "<< (toytime_stop - toytime_start)/1000000. <<" secs "<<endl;; 
				fflush(stdout);
			}
		}

		switch (test_statistics){
			case 1:
			case 2:
			case 3:
			case 31:
			case 4:
				vdata_global =  _model->GetToyData_H0();
				if(_model->hasParametricShape()){
					_model->SetTmpDataForUnbinned(_model->Get_v_pdfs_roodataset_toy());
					if(bWriteToys){
						for(int ii=0; ii<_model->Get_v_pdfs_roodataset_toy().size(); ii++){
							TString stmp = _model->Get_v_pdfs_roodataset_toy()[ii]->GetName(); stmp+="_"; stmp+=i;
							//_model->GetWorkSpace()->import(*(_model->Get_v_pdfs_roodataset_toy()[ii]), RooFit::Rename(stmp.Data()));
							wtmp->import(*(_model->Get_v_pdfs_roodataset_toy()[ii]), RooFit::Rename(stmp.Data()));
						}
					}
				}
				break;
			case 5:
			case 6:
				vdata_global = (VDChannel)_model->GetToyData_H0(_model->Get_fittedParsInData_b());
				if(_model->hasParametricShape()){
					_model->SetTmpDataForUnbinned(_model->Get_v_pdfs_roodataset_toy());
					if(bWriteToys){
						for(int ii=0; ii<_model->Get_v_pdfs_roodataset_toy().size(); ii++){
							TString stmp = _model->Get_v_pdfs_roodataset_toy()[ii]->GetName(); stmp+="_"; stmp+=i;
							//_model->GetWorkSpace()->import(*(_model->Get_v_pdfs_roodataset_toy()[ii]), RooFit::Rename(stmp.Data()));
							wtmp->import(*(_model->Get_v_pdfs_roodataset_toy()[ii]), RooFit::Rename(stmp.Data()));
						}
					}
				}
				if(!_model->UseBestEstimateToCalcQ()){
					//generate nuisance around b^hat_0 
				if(_debug>=1000)cout<<" ---- DELETEME in BuildM2lnQ_b "<<" before fluctuatenumber "<<endl;
					VChannelVSample vv =  _model->FluctuatedNumbers(0, false, 2); // toss nuisance around fitted b_hat_0 in data  
				if(_debug>=1000)cout<<" ---- DELETEME in BuildM2lnQ_b "<<" after fluctuatenumber "<<endl;
					_inputNuisances = _model->Get_randomizedPars();	
					_startNuisances = _model->Get_fittedParsInData_b();
					for(int itmp=0; itmp<_model->Get_v_pdftype().size(); itmp++){
						if(_model->Get_v_pdftype()[itmp]==typeGamma) _inputNuisances[itmp]+=1;
					}
				//	for(int i=1; i<=_model->Get_max_uncorrelation(); i++) {
				//		cout<<"DELETEMEb par "<<i<<" "<<_inputNuisances[i]<<endl;
				//	}
				}

				break;
			default:
				break;

		}
		int checkFailure = (_debug>=10?1:0);
		bool success = true;
		if(bWriteToys) Q_b[i]=0;
		else Q_b[i] = M2lnQ(success, checkFailure);
		if(success==false) {
			// skip this toy and regenerate it   --> any bias ? 
			// caveat: it may go to infinite loop if all toys fails
			i-=1;
			cout<<"WARNING: skip a failed toy, regenerating "<<endl;
		}
	}

	if(bWriteToys){
		TString stmp = "PseudoData_b_seed"; stmp+=_model->GetRdm()->GetSeed(); stmp+=".root";
		TFile *f = new TFile(stmp, "RECREATE");
		//f->WriteTObject(_model->GetWorkSpace());
		f->WriteTObject(wtmp);
		f->Close();
		return true;
	}

	if(_debug) { start_time=cur_time; cur_time=clock(); 
		cout << "\t\t\t TIME in RunMCExps run_"<<_nexps<<"_pseudo exps for b-only hypothesis: " << (cur_time - start_time)/1000. << " millisec\n";
	}

	Sort(_nexps, Q_b, iq_b, 0); // rank from small to large

	return true;
}
	void CLsBase::printM2LnQInfo(int sbANDb_bOnly_sbOnly){
		if (sbANDb_bOnly_sbOnly==0) 
			if( ( Q_b_data < Q_b[iq_b[0]] || Q_b_data > Q_sb[iq_sb[_nexps-1]] ) && _debug ){ 
				cout<<"\t probability of this -2lnQ_data is very very small, it's out of "<<_nexps<<" exps, you need more toys"<<endl;
				if(test_statistics==1){
					cout<<"\t -2lnQ_data =   "<<-2*Q_b_data+2*_nsig<<endl;
					cout<<"\t -2lnQ_b    = [ "<<-2*Q_b[iq_b[_nexps-1]]+2*_nsig<<" , "<<-2*Q_b[iq_b[0]]+2*_nsig<<" ]"<<endl;
					cout<<"\t -2lnQ_sb   = [ "<<-2*Q_sb[iq_sb[_nexps-1]]+2*_nsig<<" , "<<-2*Q_sb[iq_sb[0]]+2*_nsig<<" ]"<<endl;
				}else{
					cout<<"\t -2lnQ_data =   "<<-2*Q_b_data<<endl;
					cout<<"\t -2lnQ_b    = [ "<<-2*Q_b[iq_b[_nexps-1]]<<" , "<<-2*Q_b[iq_b[0]]<<" ]"<<endl;
					cout<<"\t -2lnQ_sb   = [ "<<-2*Q_sb[iq_sb[_nexps-1]]<<" , "<<-2*Q_sb[iq_sb[0]]<<" ]"<<endl;
				}
			}
		if(_debug>=1 && sbANDb_bOnly_sbOnly==0 ) {
			if(test_statistics==1){
				cout<<"\t -2lnQ_data =   "<<-2*Q_b_data+2*_nsig<<endl;
				cout<<"\t -2lnQ_b    = [ "<<-2*Q_b[iq_b[_nexps-1]]+2*_nsig<<" , "<<-2*Q_b[iq_b[0]]+2*_nsig<<" ]"<<endl;
				cout<<"\t -2lnQ_sb   = [ "<<-2*Q_sb[iq_sb[_nexps-1]]+2*_nsig<<" , "<<-2*Q_sb[iq_sb[0]]+2*_nsig<<" ]"<<endl;
			}else{
				cout<<"\t -2lnQ_data =   "<<-2*Q_b_data<<endl;
				cout<<"\t -2lnQ_b    = [ "<<-2*Q_b[iq_b[_nexps-1]]<<" , "<<-2*Q_b[iq_b[0]]<<" ]"<<endl;
				cout<<"\t -2lnQ_sb   = [ "<<-2*Q_sb[iq_sb[_nexps-1]]<<" , "<<-2*Q_sb[iq_sb[0]]<<" ]"<<endl;
			}
		}
		if(_debug>=100 && sbANDb_bOnly_sbOnly==0 ){
			cout<<"\t\t CHECKING order ---index Q_b Q_sb-- "<<endl;
			for(int i=0; i<_nexps; i++){
				if(test_statistics==1)	cout<<"\t "<<i<<"     "<<-2*Q_b[iq_b[i]]+2*_nsig<<"  "<<-2*Q_sb[iq_sb[i]]+2*_nsig<<endl;
				else 	cout<<"\t "<<i<<"     "<<-2*Q_b[iq_b[i]]<<"  "<<-2*Q_sb[iq_sb[i]]<<endl;
			}
		}
	}

void CLsBase::checkFittedParsInData(bool bReadPars, bool bWritePars, TString sfilename){
	_inputNuisances = _model->Get_norminalPars();	
	_startNuisances = _model->Get_norminalPars();	
	TH1D *hParsB, *hParsSB; TFile *fToWrite;
	if(bReadPars and bWritePars) { cout<<"ERROR: checkFittedParsInData():  bReadPars and bWritePars can't be both true "<<endl; exit(1); }
	if(bReadPars){
		hParsB = (TH1D*) GetTObject(sfilename.Data(),"hParsB");
		hParsSB = (TH1D*) GetTObject(sfilename.Data(),"hParsSB");
	}
	if(bWritePars){ fToWrite = new TFile(sfilename, "RECREATE");}
	// need move to a dedicated place 

	if(_model->GetTossToyConvention()==1){
		if( _model->Get_fittedParsInData_b() == 0 || bWritePars) {
			double *pars1 =new double[_model->Get_max_uncorrelation()+1];// FIXME potential memory leak
			if(bReadPars){
				for(int i=0; i<_model->Get_max_uncorrelation()+1; i++) pars1[i]=hParsB->GetBinContent(i+1); // pars1[0]: mu,  histogram starts with bin=1
				if(pars1[0]!=0) {cout<<"ERROR: checkFittedParsInData(): read hParsB, mu!=0"<<endl; exit(1);}
				if(hParsB->GetNbinsX()!=_model->Get_max_uncorrelation()+1) {
					cout<<"ERROR: checkFittedParsInData(): read hParsB->Nbins="<<hParsB->GetNbinsX()<<" != "<<
									       " model->Npars="<<_model->Get_max_uncorrelation()+1<<endl; exit(1);}
				_model->Set_fittedParsInData_b(pars1);
			}
			else {
				DoAfit(0, _model->Get_v_data_real(), _model->Get_v_pdfs_roodataset_real(), pars1);
				_model->Set_fittedParsInData_b(pars1);
				if(bWritePars){
					hParsB = new TH1D("hParsB", "hParsB", _model->Get_max_uncorrelation()+1, 0, _model->Get_max_uncorrelation()+1);
					for(int i=0; i<_model->Get_max_uncorrelation()+1; i++) hParsB->SetBinContent(i+1, pars1[i]); // pars1[0]: mu,  histogram starts with bin=1
					fToWrite->WriteTObject(hParsB);
				}
			}
		}else{
			if(bWritePars){
				hParsB = new TH1D("hParsB", "hParsB", _model->Get_max_uncorrelation()+1, 0, _model->Get_max_uncorrelation()+1);
				for(int i=0; i<_model->Get_max_uncorrelation()+1; i++) hParsB->SetBinContent(i+1, _model->Get_fittedParsInData_b()[i]); // pars1[0]: mu,  histogram starts with bin=1
				fToWrite->WriteTObject(hParsB);
			}
		}

		if( _model->Get_fittedParsInData_sb() == 0 || bWritePars) {
			double *pars2 =new double[_model->Get_max_uncorrelation()+1]; // FIXME potential memory leak
			if(bReadPars){
				for(int i=0; i<_model->Get_max_uncorrelation()+1; i++) pars2[i]=hParsSB->GetBinContent(i+1); // pars1[0]: mu,  histogram starts with bin=1
				if(pars2[0]!=_model->GetSignalScaleFactor()) {cout<<"ERROR: checkFittedParsInData(): read hParsSB, mu="<<pars2[0]<<"!= signal strength being tested: "<<_model->GetSignalScaleFactor()<<endl; exit(1);}
				if(hParsSB->GetNbinsX()!=_model->Get_max_uncorrelation()+1) {
					cout<<"ERROR: checkFittedParsInData(): read hParsSB->Nbins="<<hParsSB->GetNbinsX()<<" != "<<
									       " model->Npars="<<_model->Get_max_uncorrelation()+1<<endl; exit(1);}
				_model->Set_fittedParsInData_sb(pars2);
			}
			else {
				DoAfit(_model->GetSignalScaleFactor(), _model->Get_v_data_real(), _model->Get_v_pdfs_roodataset_real(), pars2);
				_model->Set_fittedParsInData_sb(pars2);
				if(!pars2) cout<<"pars2=0"<<endl;
				cout<<"r ="<<pars2[0]<<endl;
				if(bWritePars){
					if(_debug) cout<<" fitted sb done, saving it "<<endl;
					hParsSB = new TH1D("hParsSB", "hParsSB", _model->Get_max_uncorrelation()+1, 0, _model->Get_max_uncorrelation()+1);
					for(int i=0; i<_model->Get_max_uncorrelation()+1; i++) hParsSB->SetBinContent(i+1, pars2[i]); // pars1[0]: mu,  histogram starts with bin=1
					fToWrite->WriteTObject(hParsSB);
					if(_debug) cout<<" fitted sb done, end saving it "<<endl;
				}
			}
		}else{
			if(bWritePars){
				hParsSB = new TH1D("hParsSB", "hParsSB", _model->Get_max_uncorrelation()+1, 0, _model->Get_max_uncorrelation()+1);
				for(int i=0; i<_model->Get_max_uncorrelation()+1; i++) hParsSB->SetBinContent(i+1, _model->Get_fittedParsInData_sb()[i]); // pars1[0]: mu,  histogram starts with bin=1
				fToWrite->WriteTObject(hParsSB);
			}
		}
		if(_debug>=10){
			cout<<"*** PreFit *******, "<<_model->Get_max_uncorrelation()+1<<" pars"<<endl;
			for(int i=0; i<=_model->Get_max_uncorrelation(); i++){
				cout<<" par "<<i<<" : "<<(_model->Get_fittedParsInData_sb())[i]<<endl;
			}
		}
	}
	if(bWritePars){
		fToWrite->Close();
		//if(hParsB) delete hParsB;
		//if(hParsSB) delete hParsSB;
	}
}


bool CLsBase::BuildM2lnQ_data(){
	if(_debug)cout<<"* BuildM2lnQ_data: start"<<endl;
	_model->SetTmpDataForUnbinned(_model->Get_v_pdfs_roodataset());

	vdata_global = _model->Get_v_data();

	/*
	   if(test_statistics==1){  // otherwise it need to specify in the arguments of EvaluateLnQ() to do evaluation for data
	   _model->SetToyForUnbinned(_model->Get_v_pdfs_roodataset());
	   }
	   */

	_inputNuisances = _model->Get_norminalPars();	
	_startNuisances = _model->Get_norminalPars();	

	int checkFailure = 1;
	bool success = true;
	Q_b_data = M2lnQ(success, checkFailure, 0); // 0 for data, 1 for toy
	if(_debug)cout<<"* BuildM2lnQ_data: end "<<Q_b_data<<endl;
	return true;
}

bool CLsBase::BuildM2lnQ_sb(int nexps, bool reUsePreviousToys, bool bWriteToys){  

	_inputNuisances = _model->Get_norminalPars();	
	_startNuisances = _model->Get_norminalPars();	

	// effort for adaptive sampling
	int oldNexps = _nexps;
	vector<double> tmpQsb;
	if(reUsePreviousToys){
		// you have to make sure in the same model with same signal scale factor ...
		if(!Q_sb) reUsePreviousToys = false;  // the previous toys are either not exist or deleted
		if(nexps<=oldNexps) reUsePreviousToys = false; // if the new total nexps required is less than previous number ... 
		if(reUsePreviousToys){
			tmpQsb.clear();
			for(int i=0; i<oldNexps; i++){
				tmpQsb.push_back(Q_sb[i]);
			}
		}
	}


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


	if(Q_sb) delete [] Q_sb;
	if(iq_sb) delete [] iq_sb;
	Q_sb=new double[_nexps];
	iq_sb = new int[_nexps];	

	checkFittedParsInData();

	int nsbi, nbi;
	int tenth = _nexps/10;
	if(tenth<1) tenth=1;
	int ntemp = _nexps*_nchannels;
	if(_debug >=10 ) cout<<"ntemp="<<ntemp<<endl;
	if(test_statistics!=1 && test_statistics!=4){
		ntemp *= _model->Get_max_uncorrelation(); // if using Q_tev or Q_atlas, then multiply by the number of nuisance parameters
		ntemp *= 100;
	}
	if( ntemp>=10000000 ) {
		cout<<"\t gonna generate "<<ntemp*2<<" poisson numbers "<<endl;
	}

	clock_t toytime_start = clock(); int ntoysFor10sec = _nexps;
	for(int i=0; i<_nexps; i++){
		if(i==1) {
			clock_t toytime_stop = clock(); int timeForOneToy = toytime_stop - toytime_start; 
			if(_debug) cout<<" time per toy = "<<timeForOneToy<<" microsec"<<endl;
			if(timeForOneToy>0) ntoysFor10sec = 10*1000000/timeForOneToy+1; 
		}
		if( ntemp>=10000000 or _debug) {
			if( (i+1)%tenth == 0 ){
				printf("... Building -2lnQ,  %4.1f \%\n", i/(double)_nexps*100);
				fflush(stdout);
			}
		}
		Q_sb[i]=0;
		if(reUsePreviousToys && i<oldNexps){
			Q_sb[i] = tmpQsb[i];
			continue;
		}

		if(_debug){
			if( (i+1)%ntoysFor10sec== 0 ){
				clock_t toytime_stop = clock();
				cout<< " from_1st toy to "<<i+1<<" toy takes "<< (toytime_stop - toytime_start)/1000000. <<" secs "<<endl;; 
				fflush(stdout);
			}
		}
		switch (test_statistics){
			case 1:
			case 2:
			case 3:
			case 31:
			case 4:
				vdata_global =  _model->GetToyData_H1();
				if(_model->hasParametricShape()){
					_model->SetTmpDataForUnbinned(_model->Get_v_pdfs_roodataset_toy());
					if(bWriteToys){
						for(int ii=0; ii<_model->Get_v_pdfs_roodataset_toy().size(); ii++){
							TString stmp = _model->Get_v_pdfs_roodataset_toy()[ii]->GetName(); stmp+="_"; stmp+=i;
							_model->GetWorkSpace()->import(*(_model->Get_v_pdfs_roodataset_toy()[ii]), RooFit::Rename(stmp.Data()));
						}
					}
				}
				break;
			case 5:
			case 6:
				vdata_global = (VDChannel)_model->GetToyData_H1(_model->Get_fittedParsInData_sb());
				if(_model->hasParametricShape()){
					_model->SetTmpDataForUnbinned(_model->Get_v_pdfs_roodataset_toy());
					if(bWriteToys){
						for(int ii=0; ii<_model->Get_v_pdfs_roodataset_toy().size(); ii++){
							TString stmp = _model->Get_v_pdfs_roodataset_toy()[ii]->GetName(); stmp+="_"; stmp+=i;
							_model->GetWorkSpace()->import(*(_model->Get_v_pdfs_roodataset_toy()[ii]), RooFit::Rename(stmp.Data()));
						}
					}
				}
				if(!_model->UseBestEstimateToCalcQ()){
					VChannelVSample vv =  _model->FluctuatedNumbers(0, true, 2); // toss nuisance around fitted b_hat_mu in data 
					_inputNuisances = _model->Get_randomizedPars();	
					_startNuisances = _model->Get_fittedParsInData_sb();
					for(int itmp=0; itmp<_model->Get_v_pdftype().size(); itmp++){
						if(_model->Get_v_pdftype()[itmp]==typeGamma) _inputNuisances[itmp]+=1;
					}
					//for(int i=1; i<=_model->Get_max_uncorrelation(); i++) {
					//	cout<<"DELETEME1 par "<<i<<" "<<_inputNuisances[i]<<endl;
					//}
				}

				break;
			default:
				break;

		}
		int checkFailure = (_debug>=10?1:0);
		bool success = true;
		if(bWriteToys) Q_sb[i]=0;
		else Q_sb[i] = M2lnQ(success, checkFailure);
		if(success==false) {
			// skip this toy and regenerate it   --> any bias ? 
			// caveat: it may go to infinite loop if all toys fails
			i-=1;
			cout<<"WARNING: skip a failed toy, regenerating "<<endl;
		}
	}

	if(bWriteToys){
		TString stmp = "PseudoData_sb_seed"; stmp+=_model->GetRdm()->GetSeed(); stmp+=".root";
		TFile *f = new TFile(stmp, "RECREATE");
		f->WriteTObject(_model->GetWorkSpace());
		f->Close();
		return true;
	}
	
	if(_debug) { start_time=cur_time; cur_time=clock(); cout << "\t\t\t TIME in RunMCExps run_"<<_nexps<<"_pseudo exps: " << (cur_time - start_time)/1000. << " millisec\n"; }

	Sort(_nexps, Q_sb, iq_sb, 0);
	return true;
}

double CLsBase::M2lnQ(bool & successful, int checkFailure, int dataOrToy){
	double q = 0;
	double tmp1, tmp2, minchi2tmp;
	if(test_statistics==1){
		// change to " flutuate rates before throw H0 and then flutuate again to throw H1 ",  was "flutuate once and throw both H0 and H1"
		for(int ch=0; ch<_nchannels; ch++){
			q += (vdata_global[ch]*_lognoverb[ch]) ;
		}
		for(int ch=0; ch<_model->Get_vv_pdfs().size(); ch++){
			q += _model->EvaluateLnQ(ch, dataOrToy);// evaluate lnQ in channel i,  on the data (0) ,   1 for toy 
		}
	}else if(test_statistics==2){
		// here Q =  2ln(L_sb/L_b),  will correct in later stage to -2lnQ
		q = MinuitFit(0, tmp1, tmp1) - MinuitFit(3, tmp1, tmp1, _model->GetSignalScaleFactor());
	}else if(test_statistics==3 || test_statistics==31){
		// here Q =  2ln(L_sb/L_b),  will correct in later stage to -2lnQ
		minchi2tmp = MinuitFit(2, tmp1, tmp2);  // MinuitFit(mode, r, err_r)
		double fitted_r = tmp1;
		if(_model->AllowNegativeSignalStrength()==false && fitted_r<0) minchi2tmp = MinuitFit(0, tmp1, tmp2);  // MinuitFit(mode, r, err_r),  want r to be >=0
		if(test_statistics==3){
			if(fitted_r>=_model->GetSignalScaleFactor()) q=0;
			else q = -(MinuitFit(3, tmp1, tmp1, _model->GetSignalScaleFactor()) - minchi2tmp);
		}
		if(test_statistics==31){
			// in Feldman Cousins paper,  it allows fitted_r > the r being tested
			q = -(MinuitFit(3, tmp1, tmp1, _model->GetSignalScaleFactor()) - minchi2tmp);
		}
		if(_debug>=100)cout<<" testStat["<<test_statistics<<"]:   q = "<<q<<" fitted_r="<<fitted_r<<" minchi2tmp="<<minchi2tmp<<" tmp1="<<tmp1<<endl;
	}else if(test_statistics==4){
		// here Q =  2ln(L_sb/L_b),  will correct in later stage to -2lnQ
		minchi2tmp = MinuitFit(4, tmp1, tmp2);  // MinuitFit(mode, r, err_r)
		double fitted_r = tmp1;
		double minchi2tmp2 = 0;
		if(_model->AllowNegativeSignalStrength()==false && fitted_r<0) minchi2tmp = MinuitFit(5, tmp1, tmp2, 0);  // MinuitFit(mode, r, err_r),  want r to be >=0
		if(fitted_r>=_model->GetSignalScaleFactor()) q=0;
		else {
			minchi2tmp2 =MinuitFit(5, tmp1, tmp1, _model->GetSignalScaleFactor()); 
			q = -(minchi2tmp2 - minchi2tmp);
		}
		if(_debug>=100)cout<<" testStat["<<test_statistics<<"]: q = "<<q<<" fitted_r="<<fitted_r<<" minchi2tmp="<<minchi2tmp<<" tmp1="<<tmp1<<" minchi2tmp2="<<minchi2tmp2<<endl;
	}else if(test_statistics==5 or test_statistics==6){ // LHC type, agreed at LHC-HCG meeting on 18.05.2011   5 for upperlimit one-side, 6 for significance 
		// here Q =  2ln(L_sb/L_b),  will correct in later stage to -2lnQ

		if(_debug>=100){
			cout<<" * data in fit: ";
			for(int i=0; i<vdata_global.size(); i++){
				cout<<vdata_global[i]<<" ";
			}
			cout<<endl;

			_model->Print();
		}
		int success[1];success[0]=0;
		int success2[1];success2[0]=0;
		double * fittedPars = new double[_model->Get_max_uncorrelation()+1];
		minchi2tmp = MinuitFit(21, tmp1, tmp2, 0, fittedPars, false, checkFailure?_debug:(_debug?1:0), success);  // MinuitFit(mode, r, err_r)
		if(success[0]!=0){
			minchi2tmp = MinuitFit(21, tmp1, tmp2, 0, _model->Get_norminalPars(), true, checkFailure?_debug:(_debug?1:0), success);  // MinuitFit(mode, r, err_r)
		}
		if(success[0]!=0 and checkFailure) { 
			cout<<"ERROR WARNING data fit failed, try to dump info:  this failure sometimes related to ROOT versions, potential bugs in TMinuit. "<<endl;
			cout<<"ERROR WARNING I had experienced that 5.28.00b gave failure while 5.26 didn't on the following data card (V2011-04-21): "<<endl;
			cout<<"\n";
			cout<<"	observation 11  13 \n";
			cout<<"	bin 1 1  2  2 \n";
			cout<<"	process 0 1  0  1 \n";
			cout<<"	rate 10   50     0   50 \n";
			cout<<"	unc lnN -  2.     -   2. \n";
			cout<<endl;

			cout<<" * data in fit: ";
			for(int i=0; i<vdata_global.size(); i++){
				cout<<vdata_global[i]<<" ";
			}
			cout<<endl;
			minchi2tmp = MinuitFit(21, tmp1, tmp2, 0, 0, false, 100);  // MinuitFit(mode, r, err_r)
		}

		double fitted_r = tmp1;
		if(_model->AllowNegativeSignalStrength()==false && fitted_r<0) minchi2tmp = MinuitFit(0, tmp1, tmp2, 0, fittedPars, true);  // MinuitFit(mode, r, err_r),  want r to be >=0
		if(test_statistics==5){ // for evaluating one-sided limit 
			if(fitted_r>=_model->GetSignalScaleFactor()) q=0;
			else {
				q = -(MinuitFit(3, tmp1, tmp1, _model->GetSignalScaleFactor(), fittedPars, true, _debug?1:0, success2) - minchi2tmp);
				if(success2[0]!=0){
					q = -(MinuitFit(3, tmp1, tmp2, _model->GetSignalScaleFactor(), _model->Get_norminalPars(), true, _debug?1:0, success2) - minchi2tmp);
				}
			}

		}else if(test_statistics==6){// for evaluating significance
			q = -(MinuitFit(3, tmp1, tmp1, 0/*fixed at mu=0*/, fittedPars, true, _debug?1:0, success2) - minchi2tmp);
			if(success2[0]!=0){
				q = -(MinuitFit(3, tmp1, tmp2, 0, _model->Get_norminalPars(), true, _debug?1:0, success2) - minchi2tmp);
			}
		}
		if(_debug>=100)cout<<" testStat["<<test_statistics<<"]: q = "<<q<<" fitted_r="<<fitted_r<<" minchi2tmp="<<minchi2tmp<<" tmp1="<<tmp1<<endl;

		if(checkFailure){
			if(fitted_r>=_model->GetSignalScaleFactor()){
				if(_debug)cout<<"data OverFlow:  fitted_r= "<<fitted_r<<",  the probe r ="<<_model->GetSignalScaleFactor()<<endl;
			}
			if(fitted_r<0){
				if(_debug)cout<<"data UnderFlow:  fitted_r= "<<fitted_r<<",  the probe r ="<<_model->GetSignalScaleFactor()<<endl;
			}
			if(_debug>=100) cout<<" end of data fit"<<endl;
		}
		if(fittedPars)delete [] fittedPars;

		if((success2[0]!=0 or success[0]!=0) && !checkFailure) cout<<"FAILED_TOY : fit unconverged"<<endl;
		if((success2[0]!=0 or success[0]!=0) && checkFailure) cout<<"FAILED_DATA : fit unconverged"<<endl;
		if((success2[0]==0 and success[0]==0) && !checkFailure && _debug) cout<<"SUCCESSFUL_TOY : fit converged"<<endl;
		if((success2[0]==0 and success[0]==0) && checkFailure && _debug) cout<<"SUCCESSFUL_DATA : fit converged"<<endl;
		if(success[0]!=0 or success2[0]!=0) successful = false;
	}
	if(_debug>=100) cout<<"-2lnQ = "<<q<<endl;
	return q;
}
void CLsBase::prepareLogNoverB(){ // only necessary when evaluating LEP type statistics 
	_nsig=0; _nbkg=0; _ndat=0;
	vector<double> vs, vb, vd; 
	vs.clear(); vb.clear(); vd.clear();
	_nchannels = _model->NumOfChannels();
	double br = _model->GetSignalScaleFactor();
	const VChannelVSample & vv = _model->Get_vv_exp_sigbkgs();
	for(int i=0; i<_nchannels; i++){
		double totbkg = 0, totsig=0;
		for(int isamp = 0; isamp<(_model->Get_vv_exp_sigbkgs())[i].size(); isamp++){
			if(_debug>=100) cout<<"ch "<<i<<" isamp "<<isamp<<",  nsigproc= "<<_model->GetNSigprocInChannel(i)<<endl;
			if(isamp<_model->GetNSigprocInChannel(i)) totsig+= (_model->Get_vv_exp_sigbkgs())[i][isamp];
			else totbkg+=(_model->Get_vv_exp_sigbkgs())[i][isamp];
		}

		if(_model->GetPhysicsModel()==typeChargedHiggs){
			// ************** for charged higgs    generalized for all channels
			double tmp = 0;
			// HH_ltau
			tmp +=  (br*br*vv[i][0]); // HH -> r^2 * HH
			// WH_ltau
			tmp +=  (2*br*(1-br)*vv[i][1]); //WH -> 2*r*(1-r)*WH
			// WW_ltau (tt->ltau, real tau) 
			tmp +=  ((1-br)*(1-br)*vv[i][2]); //WW -> (1-r)*(1-r)*WW
			// WW_ltau (tt~->ll, also part of WW, one lepton fakes as tau) 
			tmp +=  ((1-br)*(1-br)*vv[i][3]); //WW -> (1-r)*(1-r)*WW

			tmp -= vv[i][2];
			tmp -= vv[i][3];

			totsig = tmp;
			// -----------------end  for charged higgs
			cout<<" BRANCH ratio = "<<br<<" totsig in channel ["<<i<<"]:  totsig = "<<totsig<<endl;
		}


		vs.push_back(totsig);
		vb.push_back(totbkg);
		vd.push_back((_model->Get_v_data())[i]);
		_nsig += vs[i];
		_nbkg += vb[i];
		_ndat += vd[i];
	}

	double n, noverb;
	if(_lognoverb) delete [] _lognoverb;
	_lognoverb=new double[_nchannels];
	for(int i=0; i<_nchannels; i++){	
		// skip a channle in which nsig==0 || ntotbkg==0
		if(( _model->AllowNegativeSignalStrength()==true || vs[i] > 0) && vb[i] > 0) {
			n=vs[i]+vb[i];
			noverb=n/vb[i];
			_lognoverb[i]= ( n>0 ?log(noverb):0 );
			if(_model->GetPhysicsModel()==typeStandardModel)_lognoverb[i]=fabs(_lognoverb[i]); // FIXME  treat downward fluctuation properly for excess/deficit signal. for SM higgs, we look for excess, while for Charged higgs, they may look for deficit, please take the fabs out
			if(_model->GetPhysicsModel()==typeChargedHiggs)_lognoverb[i]=_lognoverb[i]; // FIXME  treat downward fluctuation properly for excess/deficit signal. for SM higgs, we look for excess, while for Charged higgs, they may look for deficit, please take the fabs out
		}else{
			_lognoverb[i]=0;
		}
		if(_debug>=10)cout<<" \t channel "<<i<<" s="<<vs[i]<<" b="<<vb[i]<<" d="<<vd[i]<<" lognoverb="<<_lognoverb[i]<<endl;	
	}
}

double CLsBase::FindLimitFromPreComputedGrid(std::map<double, TTree*> gridCLsb, std::map<double, TTree*> gridCLb, std::map<double, double> gridQdata, double alpha, TString plotName){ // from precomputed m2lnQ grid to extract r corresponding to _alpha ... e.g. 0.05
	TGraphErrors *tge =  new TGraphErrors();
	if (_debug >= 10) std::cout << "Search for upper limit using pre-computed grid of p-values" << std::endl;

	int i = 0, n = gridCLsb.size();
	if(_debug)cout<<" grid size = "<<n<<endl;
	for (std::map<double, TTree *>::iterator itg = gridCLsb.begin(), edg = gridCLsb.end(); itg != edg; ++itg) {
		double cls, clserr;
		GetCLs(gridQdata[itg->first], itg->second, gridCLb[itg->first], cls, clserr, _debug);
		tge->SetPoint(     i, itg->first, cls   ); 
		tge->SetPointError(i, 0,          clserr);
		if(_debug>=10)cout<<" input grid:  r="<<itg->first<<" cls="<<cls<<"+/-"<<clserr<<endl;
		i++;
	}

	double limit, limitErr;
	FindLimitFromTGE(tge, alpha, limit, limitErr, plotName);

	delete tge;
	return limit;
}

double CLsBase::FindLimitFromTGE(TGraphErrors *tge, double alpha, double &limit, double & limitErr, TString plotName){
	tge->Sort();
	TF1 *fit = fitToRvsCL_expo;

	double dist=9999.;
	int n= tge->GetN();
	double rMin=0, rMax=0;
	std::pair<double, double> clsMin, clsMax;
	for (int i = 0; i < n; i++) {
		double x = tge->GetX()[i], y = tge->GetY()[i], ey = tge->GetErrorY(i);
		if (y-3*ey >= alpha) { rMin = x; clsMin = make_pair(y,ey); }
		if (y+3*ey <= alpha) { rMax = x; clsMax = make_pair(y,ey); }
		if (fabs(y-alpha) < dist) { limit = x; dist = fabs(y-alpha); }
		if (_debug >=10 ) std::cout << "  r " << x << ", CLs = " << y << " +/- " << ey << std::endl;
	}
	if((clsMin.first==0 and clsMin.second==0) || (clsMax.first==0 and clsMax.second==0))
	{
		if(_debug) cout<<"Couldn't find both points with CLs < alpha-3*eCLs and CLs > alpha+3*eCLs, decide to fit to full available range"<<endl;
		rMin = tge->GetX()[0]; clsMin=make_pair(tge->GetY()[0], tge->GetErrorY(0)); 
		rMax = tge->GetX()[n-1]; clsMax=make_pair(tge->GetY()[n-1], tge->GetErrorY(n-1));
	}

	limitErr = std::max(limit-rMin, rMax-limit);
	fit->SetRange(rMin,rMax);

	if (_debug) {
		std::cout << " Limit Before Fit: r = " << limit << " +/- " << limitErr <<endl;
		std::cout << " r fitting range: [" << rMin << ", " << rMax << "]"<<endl;
	}

	fit->FixParameter(0,alpha);
	fit->SetParameter(1,log(clsMax.first/clsMin.first)/(rMax-rMin));
	fit->SetParameter(2,limit);
	double rMinBound, rMaxBound; fit->GetRange(rMinBound, rMaxBound);
	limitErr = std::max(fabs(rMinBound-limit), fabs(rMaxBound-limit));
	int npoints = 0; 
	for (int j = 0; j < tge->GetN(); ++j) { 
		if (tge->GetX()[j] >= rMinBound && tge->GetX()[j] <= rMaxBound) npoints++; 
	}

	tge->Fit(fit,(_debug <= 10 ? "QNR EX0" : "NR EXO"));
	limit = fit->GetParameter(2);
	limitErr = fit->GetParError(2);
	if (_debug) {
		std::cout << "After Fit to " << npoints << " points: r@"<<(1-alpha)*100<<"%CL = " << limit << " +/- " << limitErr << std::endl;
	}

	if (plotName!="") {
		//FIXME to-dos: change background to white,  add x/y labels,  add markers for central values .... show limit+/-err on the plot 
		TCanvas c1("c1","c1");
		tge->Sort();
		tge->SetLineWidth(2);
		tge->Draw("AP");
		fit->Draw("SAME");
		TLine line(tge->GetX()[0], alpha, tge->GetX()[tge->GetN()-1], alpha);
		line.SetLineColor(kRed); line.SetLineWidth(2); line.Draw();
		line.DrawLine(limit, 0, limit, tge->GetY()[0]);
		line.SetLineWidth(1); line.SetLineStyle(2);
		line.DrawLine(limit-limitErr, 0, limit-limitErr, tge->GetY()[0]);
		line.DrawLine(limit+limitErr, 0, limit+limitErr, tge->GetY()[0]);
		c1.Print(plotName+".gif");
		c1.Print(plotName+".root");
	}

	return limit;
}
};
