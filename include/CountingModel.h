#ifndef CountingModel_H
#define CountingModel_H
/*
 * =====================================================================================
 *
 *       Filename: CountingModel.h 
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  03/10/2010 10:50:40 PM CET
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Mingshui Chen (), Mingshui.Chen@cern.ch
 *        Company:  Univsity of Florida
 *
 * =====================================================================================
 */
#include <vector>
#include <string>
#include "CRandom.h"
#include "Utilities.h"
#include "SMHiggsBuilder.h"
#include "RooAbsPdf.h"
#include "RooRealVar.h"
#include "RooBifurGauss.h"
#include "RooWorkspace.h"
#include <map>
#include "TString.h"
#include "TH1D.h"

using namespace std;
namespace lands{
	typedef vector<double> VDChannel;
	//typedef vector<int> VIChannel;
	typedef vector<double> VIChannel;
	typedef vector< vector<double> > VChannelVSample;
	typedef vector< vector< vector<int> > > VChannelVSampleVUncertainty;
	typedef vector< vector< vector< vector<double> > > > VChannelVSampleVUncertaintyVParameter;
	typedef vector< vector< vector< vector<float> > > > VChannelVSampleVsetVval;
	typedef map< string, vector< vector<double> > > MapStrVV;
	typedef map< string, vector<double> > MapStrV;

	typedef vector<TH1D*> VChannelTH;
	typedef vector< vector<TH1D*> > VChannelVSampleTH;
	typedef vector< vector< vector< vector<TH1D*> > > > VChannelVSampleVUncertaintyVParameterTH;
	typedef map< string, vector< vector<TH1D*> > > MapStrVVTH;
	typedef map< string, vector<TH1D*> > MapStrVTH;

	enum enumPdfType {typeLogNormal=1, typeTruncatedGaussian=2, typeGamma=3, typeShapeGaussianLinearMorph=4, typeShapeGaussianQuadraticMorph=5, 
		typeBifurcatedGaussian=6, typeFlat=7, typeControlSampleInferredLogNormal=11 };
	// for uncertainties only affecting Shape, we can choose different morphing algorithm.  in commom, user must provide three templates:  norminal,  shift_1sigma_up, shift_1sigma_down
	// i.e. 3 parameters for each shape uncertainty in each bin .... ,  the interface alway do normalization (to unity) for all three templates.
	
	enum enumPhysicsModel {typeModelBegin=0,typeStandardModel=1, typeChargedHiggs=2, typeCvCfHiggs=3, typeC5Higgs=4, typeModelEnd};

	enum enumDecayMode {decayHZZ=11, decayHWW=10, decayHTT=2, decayHBB=1, decayHGG=8, decayHZG=9, decayHCC=5, decayHTopTop=6, decayHGluGlu=7, decayHSS=4, decayHMM=3}; // add more later
	enum enumProductionMode {productionGGH=1, productionVH=2, productionQQH=3, productionTTH=4}; // VH contains WH and ZH,    QQH aka VBF
	struct structPOI {
		TString name;
		double value;
		double errUp;
		double errDown;
		double minV;
		double maxV;
		structPOI(TString s, double v, double eu, double ed, double minv, double maxv):name(s),value(v),errUp(eu),errDown(ed),minV(minv), maxV(maxv){}
	};

	class CountingModel
	{ 
		public:
			CountingModel();
			~CountingModel();
			void AddChannel(std::string channel_name, double num_expected_signal, double num_expected_bkg_1, double num_expected_bkg_2=-1, 
						double num_expected_bkg_3=-1, double num_expected_bkg_4=-1, double num_expected_bkg_5=-1, double num_expected_bkg_6 = -1 );
			void AddChannel(double num_expected_signal, double num_expected_bkg_1, double num_expected_bkg_2=-1, 
						double num_expected_bkg_3=-1, double num_expected_bkg_4=-1, double num_expected_bkg_5=-1, double num_expected_bkg_6 = -1 );

			void AddChannel(double num_expected_signal, vector<double> num_expected_bkgs);
			void AddChannel(std::string channel_name, double num_expected_signal, vector<double> num_expected_bkgs);
			void AddChannel(vector<double> num_expected_yields, int signal_processes=1);
			void AddChannel(std::string channel_name, vector<double> num_expected_yields, int signal_processes=1, int decaymode=-1);
			void AddChannel(std::string channel_name, vector<double> num_expected_signals, vector<double> num_expected_bkgs, int decaymode=-1);

			void TagUncertaintyFloatInFit(string uncname, bool b);
			void TagUncertaintyFloatInFit(int uncIndex, bool b);

			// LogNormal and TruncatedGaussian 
			void AddUncertainty(int index_channel, int index_sample, double uncertainty_in_relative_fraction, int pdf_type, std::string uncname );
			void AddUncertainty(string chname, int index_sample, double uncertainty_in_relative_fraction, int pdf_type, std::string uncname );
			void AddUncertainty(int index_channel, int index_sample, double uncertainty_in_relative_fraction, int pdf_type, int index_correlation );
			// for asymetric uncertainties 
			void AddUncertainty(int index_channel, int index_sample, double uncertainty_in_relative_fraction_down, double uncertainty_in_relative_fraction_up, int pdf_type, std::string uncname );
			void AddUncertainty(string chname, int index_sample, double uncertainty_in_relative_fraction_down, double uncertainty_in_relative_fraction_up, int pdf_type, std::string uncname );
			void AddUncertainty(int index_channel, int index_sample, double uncertainty_in_relative_fraction_down, double uncertainty_in_relative_fraction_up, int pdf_type, int index_correlation);

			// From SideBand
			void AddUncertainty(int index_channel, int index_sample, double rho, double rho_err, double B, int pdf_type, std::string uncname );
			void AddUncertainty(string chname, int index_sample, double rho, double rho_err, double B, int pdf_type, std::string uncname );
			void AddUncertainty(int index_channel, int index_sample, double rho, double rho_err, double B, int pdf_type, int index_correlation );
			// when B is large, it's close to Gaussian. then use Gaussian, don't use Gamma function

			// For Shape parameters
			void AddUncertainty(int index_channel, int index_sample, int npar, double *par, int pdf_type, int index_correlation );
			void AddUncertainty(int index_channel, int index_sample, int npar, double *par, int pdf_type, std::string uncname );
			void AddUncertainty(string chname, int index_sample, int npar, double *par, int pdf_type, std::string uncname );

			void AddObservedData(int index_channel, double num_data);
			void AddObservedData(string c, double num_data);
			void SetData(VDChannel data, bool bRealData=true){v_data=data; if(bRealData)v_data_real=data;}
			void SetData(TString filename, TString datahistname, bool bRealData=true);
			void SetData(vector<int> data){for(int i=0; i<data.size(); i++) v_data[i]=data[i];}

			const VChannelVSample& Get_vv_exp_sigbkgs_nonscaled() {return vv_exp_sigbkgs;}
			const VChannelVSample& Get_vv_exp_sigbkgs()           {return vv_exp_sigbkgs_scaled;}
			void Set_vv_randomized_sigbkgs(const VChannelVSample& vv){vv_randomized_sigbkgs=vv;}
			const VChannelVSample& Get_vv_randomized_sigbkgs(){return vv_randomized_sigbkgs_scaled;}
			void Set_vv_fitted_sigbkgs(const VChannelVSample& vv){vv_fitted_sigbkgs=vv;}
			const VChannelVSample& Get_vv_fitted_sigbkgs(){return vv_fitted_sigbkgs;}
			void Set_vv_fitted_sigbkgs_scaled(const VChannelVSample& vv){vv_fitted_sigbkgs_scaled=vv;}
			const VChannelVSample& Get_vv_fitted_sigbkgs_scaled(){return vv_fitted_sigbkgs_scaled;}
			const VDChannel& Get_v_data(){return v_data;}
			VDChannel Get_AsimovData(int b);
			const VDChannel& Get_v_data_real(){return v_data_real;}
			const VDChannel& Get_v_exp_sigbkgs(int channel){return vv_exp_sigbkgs_scaled[channel];}

			double GetExpectedNumber(int index_channel, int index_sample);

			void Print(int printLevel=0);
			bool Check();

			VChannelVSample FluctuatedNumbers(double *pars = NULL, bool scaled=true, int bUseBestEstimateToCalcQ=1, bool includeCountingParts=true); 
			// 0 use randomized set,  1 use best estimates a priori, 2 use fitted posterior (after looking at data)
			
			VIChannel GetToyData_H0(double *pars=NULL);// background only hypothesis
			VIChannel GetToyData_H1(double *pars=NULL);// alternative hypothesis
			
			int NumOfChannels(){return vv_exp_sigbkgs.size();}
			void SetUseSystematicErrors(bool b){b_systematics=b;ConfigUncertaintyPdfs();}
			bool IsUsingSystematicsErrors(){return b_systematics;}

			void SetSignalScaleFactor(double r, int bScaleBestEstimate=1);
			double GetSignalScaleFactor(){return _common_signal_strength; }
			void SetRdm(CRandom *rdm){_rdm=rdm;}
			CRandom * GetRdm(){return _rdm;}

			void UseAsimovData(int b=0);


			const VChannelVSampleVUncertaintyVParameter& Get_vvvv_uncpar(){return vvvv_uncpar;}
			void Set_vvvv_uncpar(const VChannelVSampleVUncertaintyVParameter& vvvv){vvvv_uncpar=vvvv;}
			const VChannelVSampleVUncertainty& Get_vvv_pdftype(){return vvv_pdftype;}
			const VChannelVSampleVUncertainty& Get_vvv_idcorrl(){return vvv_idcorrl;}

			int RemoveChannelsWithExpectedSignal0orBkg0(int king = 2); // king=0 for removing a bin with bkg=0, 1 for sig=0, 2 for bkg=0||sig=0  

			int Get_max_uncorrelation() {return max_uncorrelation;}
			const VDChannel& Get_v_TruncatedGaussian_maxUnc() {return v_TruncatedGaussian_maxUnc;}
			const vector<int>& Get_v_pdftype() {return v_pdftype;}

			void SetAllowNegativeSignalStrength(bool b){b_AllowNegativeSignalStrength = b;}
			bool AllowNegativeSignalStrength(){return b_AllowNegativeSignalStrength;}

			const vector<double>& Get_v_GammaN(){return v_GammaN;}
			const vector<std::string>& Get_v_uncname(){return v_uncname;}
			const vector<bool>& Get_v_uncFloatInFit(){return v_uncFloatInFit;}
			void SetDebug(int i){_debug=i;}
			int GetDebug(){return _debug;}
			std::string GetChannelName(int i){return v_channelname.at(i);} //need make check 
			int GetNSigprocInChannel(int i){return v_sigproc.at(i);} //need make check
			const vector<int>& Get_v_sigproc(){return v_sigproc;}

			MapStrV Get_map_flatPars(){return map_flatPars;}
			void Set_flatPars(pair<string, vector<double> > f);
			void SetFlatParInitVal(string s, double d);

			void SetProcessNames(int ch, vector<std::string> vproc); //need make check
			void SetProcessNames(string ch, vector<std::string> vproc); //need make check
			const vector<std::string>& GetProcessNames(int ch){return vv_procname[ch];}  //need make check
			void SetModelName(const std::string& s);
			const std::string& GetModelName(){return _modelName;}
			int GetDecayMode(){return _decayMode;}
			void MakeListOfShapeUncertainties();
			const vector<int>& GetListOfShapeUncertainties(int c, int p){ return vvv_shapeuncindex[c][p]; } // need make check
			void SetMoveUpShapeUncertainties(bool b){bMoveUpShapeUncertainties=b;}
			bool GetMoveUpShapeUncertainties(){return bMoveUpShapeUncertainties;}

			void SetTossToyConvention(int c){_tossToyConvention = c;}
			int GetTossToyConvention(){return _tossToyConvention;}
			double* Get_fittedParsInData_sb(){return _fittedParsInData_sb;}
			double* Get_fittedParsInData_b(){return _fittedParsInData_bonly;}
			void Set_fittedParsInData_sb(double *p){
			       	_fittedParsInData_sb=p;
				vv_fitted_sigbkgs_scaled = FluctuatedNumbers(p); // fit_sb and scaled
				vv_pdfs_norm_fitted_scaled = vv_pdfs_norm_varied;
			}
			void Set_fittedParsInData_b(double *p){
				_fittedParsInData_bonly=p;
				vv_fitted_sigbkgs= FluctuatedNumbers(p, 0);
				vv_pdfs_norm_fitted = vv_pdfs_norm_varied;
			}

			void SetUseBestEstimateToCalcQ(int b=1){ _UseBestEstimateToCalcQ=b;}
			int UseBestEstimateToCalcQ(){ return _UseBestEstimateToCalcQ;}

			// start to add parametric shape into model
			void SetObservableRange(int n, double xmin, double xmax){bRedefineObservableRange = true; ObservableBins=n; ObservableRangeMin=xmin; ObservableRangeMax=xmax;};

			bool hasParametricShape(){return bHasParametricShape;}
			const vector< vector<string> >& Get_vv_pdfs(){return vv_pdfs;}
			const vector< vector<double *> >& Get_vv_pdfs_params(){return vv_pdfs_params;} //nominal parameters for each pdf
			const vector< vector<double> >& Get_vv_pdfs_norm_nonscaled(){return vv_pdfs_norm;} //normalization of each pdf, i.e. expected number of events in each process 
			const vector< vector<int> >& Get_vv_pdfs_npar(){return vv_pdfs_npar;} // number of parameters for each pdf
			const vector< int >& Get_v_pdfs_nbin(){return v_pdfs_nbin;} // need for throwing random numbers,  should consistent for all pdfs in a channel, by default 100 
			const vector< double >& Get_v_pdfs_xmin(){return v_pdfs_xmin;}
			const vector< double >& Get_v_pdfs_xmax(){return v_pdfs_xmax;}
			const vector< vector<double> >& Get_vv_pdfs_data(){return vv_pdfs_data;} // in each channel, it has a list of events
			// uncertainties ....   each source affects parameters  -->  two additional sets of parameters
			const vector< vector< vector<int> > >& Get_vvv_pdfs_idcorrl(){return vvv_pdfs_idcorrl;}
			const vector< vector< vector<int> > >& Get_vvv_pdfs_pdftype(){return vvv_pdfs_pdftype;}
			// three types of uncertainties: 1. only affect shape;  2. only affect normalization; 3. affect both 
			const vector< vector< vector<int> > >& Get_vvv_pdfs_unctype(){return vvv_pdfs_unctype;}
			const vector< vector< vector<double*> > >& Get_vvv_pdfs_params_up(){return vvv_pdfs_params_up;} // correponding to up shift of the uncertainty source, e.g. jet energy scale
			const vector< vector< vector<double*> > >& Get_vvv_pdfs_params_down(){return vvv_pdfs_params_down;} //  correponding to down shift of the uncertainty source, e.g. jet energy scale
			const vector< vector< vector< vector<double> > > >& Get_vvv_pdfs_normvariation(){return vvv_pdfs_normvariation;} // correponding to normalization changes of the uncertainty source
			const vector< vector<double *> >& Get_vv_pdfs_params_varied(){return vv_pdfs_params_varied;} //nominal parameters for each pdf
			const vector< vector<double> >& Get_vv_pdfs_norm_varied(){return vv_pdfs_norm_varied;} //normalization of each pdf, i.e. expected number of events in each process 
			const vector< vector<double> >& Get_vv_pdfs_norm_scaled(){return vv_pdfs_norm_scaled;} //normalization of each pdf, i.e. expected number of events in each process 
			const vector< vector<double> >& Get_vv_pdfs_data_toy(){return vv_pdfs_data_toy;} // in each channel, it has a list of events
			void Set_vv_pdfs_norm_randomized(const VChannelVSample& vv){vv_pdfs_norm_randomized=vv;}
			const VChannelVSample& Get_vv_pdfs_norm_randomized(){return vv_pdfs_norm_randomized_scaled;}
			void Set_vv_pdfs_norm_fitted(const VChannelVSample& vv){vv_pdfs_norm_fitted=vv;}
			const VChannelVSample& Get_vv_pdfs_norm_fitted(){return vv_pdfs_norm_fitted;}
			void Set_vv_pdfs_norm_fitted_scaled(const VChannelVSample& vv){vv_pdfs_norm_fitted_scaled=vv;}
			const VChannelVSample& Get_vv_pdfs_norm_fitted_scaled(){return vv_pdfs_norm_fitted_scaled;}
			const vector<int>& Get_v_pdfs_sigproc(){return v_pdfs_sigproc;}
			const vector< RooAbsData* >& Get_v_pdfs_roodataset();
			const vector< RooAbsData* >& Get_v_pdfs_roodataset_toy(){return v_pdfs_roodataset_toy;} // in each channel, it has a list of events
			const vector< RooAbsData* >& Get_v_pdfs_roodataset_real(){return v_pdfs_roodataset_real;} // in each channel, it has a list of events
			const vector< RooAbsData* >& Get_v_pdfs_roodataset_asimovb(){return v_pdfs_roodataset_asimovb;} // in each channel, it has a list of events
			void ConstructAsimovData(int b, bool nominal=true, double injectMu =1.);
			const vector< RooAbsData* >& Get_v_pdfs_roodataset_asimovsb(){return v_pdfs_roodataset_asimovsb;} // in each channel, it has a list of events
			const vector< double >& Get_v_pdfs_floatParamsVaried(){return v_pdfs_floatParamsVaried;}
			const vector< vector< double > >& Get_v_pdfs_floatParamsUnc(){ return v_pdfs_floatParamsUnc;}
			void Set_v_pdfs_floatParamsUnc(const vector< vector< double > >& vv ){ v_pdfs_floatParamsUnc = vv;}
			void SetFlatParameterRange(int id, double middle, double errLow, double errUp);
			void SetFlatNormalizationRange(int id, double errLow, double errUp);
			const vector<int>& Get_v_pdfs_floatParamsIndcorr() {return v_pdfs_floatParamsIndcorr;}      // only for params
			const vector<int>& Get_v_pdfs_floatParamsType() {return v_pdfs_floatParamsType;}      // only for params

			const vector<string>& Get_v_pdfs_channelname(){return v_pdfs_channelname;}
			const vector< vector<string> >& Get_vv_pdfs_procname(){return vv_pdfs_procname;}
			const vector<TString>& Get_v_pdfs_sb(){return v_pdfs_sb;}
			const vector<TString>& Get_v_pdfs_s(){return v_pdfs_s;}
			const vector<TString>& Get_v_pdfs_b(){return v_pdfs_b;}
			const vector<TString>& Get_v_pdfs_observables(){return v_pdfs_observables;}
			RooWorkspace * GetWorkSpace(){return _w;}
			RooWorkspace * GetWorkSpaceVaried(){return _w_varied;}

			const vector< vector< vector<int> > >& Get_vvv_pdfs_nuisancesindex(){return vvv_pdfs_nuisancesindex;}
			const VChannelVSampleVUncertaintyVParameter& Get_vvvv_pdfs_ChProcSetEvtVals(){return vvvv_pdfs_ChProcSetEvtVals;}
			void Set_vvvv_pdfs_ChProcSetEvtVals(const VChannelVSampleVUncertaintyVParameter& vvvv_pdfs){vvvv_pdfs_ChProcSetEvtVals=vvvv_pdfs;}
			const VChannelVSampleVUncertaintyVParameter& Get_vvvv_pdfs_ChProcSetParVals(){return vvvv_pdfs_ChProcSetParVals;}
			void Set_vvvv_pdfs_ChProcSetParVals(const VChannelVSampleVUncertaintyVParameter& vvvv_pdfs){vvvv_pdfs_ChProcSetParVals=vvvv_pdfs;}

			void AddChannel(string channel_name, vector<string> vprocname , RooRealVar* observable, vector<RooAbsPdf*> sigPdfs, vector<double> sigNorms, vector<RooAbsArg*> vsExtraNorm, 
				       	vector<RooAbsPdf*> bkgPdfs, vector<double> bkgNorms, vector<RooAbsArg*> vbExtraNorm, int decaymode=-1);
			// need to add names of each parameter .... 
			double EvaluateLnQ(int ch, int dataOrToy); // for Likelihood ratio
			double EvaluateChi2(double *par, vector<float>& v_cachPdfValuestmp, vector< vector< vector<float> > >& vvv_cachPdfValuestmp, int bUseBestEstimateToCalcQ=1);          // for Chi2
			double EvaluateGL(int ch, double xr); // for bayesian 
			double EvaluateGL(vector< vector<double> > vvnorms, vector<double> vparams, double xr, VChannelVSample& vvs, VChannelVSample&vvb); // for bayesian 
			void AddObservedDataSet(int index_channel, RooAbsData* rds);
			void AddObservedDataSet(string channelname, RooAbsData* rds);
			void SetDataForUnbinned(vector< RooAbsData*> data, bool bRealData=true);
			void SetDataForUnbinned(TString filename, bool bRealData=true);
			void SetTmpDataForUnbinned(vector< RooAbsData*> data); // set it before doing fit
			void SetToyForUnbinned(vector< RooAbsData*> data); // set it before doing fit
			void AddUncertaintyOnShapeNorm(int index_channel, int index_sample, double uncertainty_in_relative_fraction_down, double uncertainty_in_relative_fraction_up, int pdf_type, int index_correlation );
			void AddUncertaintyOnShapeNorm(int index_channel, int index_sample, double uncertainty_in_relative_fraction_down, double uncertainty_in_relative_fraction_up, int pdf_type, string uncname);
			void AddUncertaintyOnShapeNorm(string chname, int index_sample, double uncertainty_in_relative_fraction_down, double uncertainty_in_relative_fraction_up, int pdf_type, string uncname);
			// From SideBand
			void AddUncertaintyOnShapeNorm(int index_channel, int index_sample, double rho, double rho_err, double B, int pdf_type, std::string uncname );
			void AddUncertaintyOnShapeNorm(string chname, int index_sample, double rho, double rho_err, double B, int pdf_type, std::string uncname );
			void AddUncertaintyOnShapeNorm(int index_channel, int index_sample, double rho, double rho_err, double B, int pdf_type, int index_correlation );

			bool AddUncertaintyOnShapeParam(string pname, double mean, double sigmaL, double sigmaR, double rangeMin=0, double rangeMax=0 );
			bool AddUncertaintyOnShapeParam(string pname);// add flatParam,  values(norminal, range) taken from workspace 
			bool AddFlatParam(string pname, double norminal, double rangeMin, double rangeMax);// add flatParam, values taken from text file

			void AddUncertaintyAffectingShapeParam(string uname, string pname, double sigmaL, double sigmaR);

			const MapStrVV& Get_map_param_sources(){return map_param_sources;}

			double * Get_norminalPars(){return _norminalPars;}
			double * Get_norminalParsSigma(){return _norminalParsSigma;}
			double * Get_randomizedPars(){return _randomizedPars;}

			void SetMass(double d);

			const vector< vector<TString> >& Get_vv_pdfs_extranormNAME() {return vv_pdfs_extranormNAME;}

			void FlagChannelsWithParamsUpdated(int i);
			void UnFlagAllChannels(bool b);
			void FlagAllChannels();
			const vector< vector<bool> > & Get_vv_statusUpdated(){return vv_statusUpdated;};

			void Set_minuitSTRATEGY(int i){minuitSTRATEGY=i;};
			int Get_minuitSTRATEGY(){return minuitSTRATEGY;};
			void Set_minuitTolerance(double i){minuitTolerance=i;};
			int Get_minuitTolerance(){return minuitTolerance;};
			void Set_nuisancesRange(double i){nuisancesRange=i;};
			int Get_nuisancesRange(){return nuisancesRange;};
			void Set_maximumFunctionCallsInAFit(int i){maximumFunctionCallsInAFit=i;};
			int Get_maximumFunctionCallsInAFit(){return maximumFunctionCallsInAFit;};

			void SetPhysicsModel(int i){
				_PhysicsModel=i;
				if(i<=typeModelBegin && i>=typeModelEnd)
				{ cout<<"ERROR: we only support typeStandardModel and typeChargedHiggs, typeCvCfHiggs, typeC5Higgs"<<endl; exit(1);}
				if(i==typeCvCfHiggs or i==typeC5Higgs){ _smhb = new SMHiggsBuilder(); }
			};
			int GetPhysicsModel(){return _PhysicsModel;};

			void ForceSymmetryError(bool b){b_ForceSymmetryError = b;};
			void MultiSigProcShareSamePDF(bool b){b_MultiSigProcShareSamePDF = b;};


			void DumpFitResults(double *pars, TString ssave);
			void Set_maxSetsCaching(int i){maxsets_forcaching=i;};
			void Set_PrintParameter(int i, int j){PrintParameterFrom =i; PrintParameterTo = j;};
			int GetPrintParameterFrom(){return PrintParameterFrom;};
			int GetPrintParameterTo(){return PrintParameterTo;};

			int Get_MH_i(){return _MH_i;};

			bool bFixingPOIs(){return _bFixingPOIs;};
			const vector< std::pair<TString, double> > & GetPOIsToBeFixed(){return vsPOIsToBeFixed;}
			void SetbFixingPOIs(bool b){_bFixingPOIs = b; if(!b) vsPOIsToBeFixed.clear();};
			void SetPOItoBeFixed(TString sp, double p){ vsPOIsToBeFixed.push_back(std::make_pair(sp,p)); _bFixingPOIs = true;};
			const vector<structPOI> & POIs(){return vPOIs;};
			void addPOI(structPOI p) {bool b=false; for(int i=0; i<vPOIs.size(); i++) {if(vPOIs[i].name==p.name) b=true;} if(!b)vPOIs.push_back(p); };
			void setPOI(int i, double v, double eu, double ed){vPOIs[i].value=v; vPOIs[i].errUp=eu; vPOIs[i].errDown=ed;};
			void keepOnlyPOI(TString spoi);
			void Set_minuitPrintLevel(int i){minuitPrintLevel = i;};
			int Get_minuitPrintLevel(){return minuitPrintLevel;};

			void SetNoErrorEstimation(bool b){noErrorEstimation=b;};
			bool GetNoErrorEstimation(){return noErrorEstimation;};

			void AddCouplingParameter(TString s);
			void CheckCouplingSet();

		const	vector< vector<int> > & Get_vv_channelDecayMode(){return vv_channelDecayMode;};
		const	vector< vector<int> > & Get_vv_pdfs_channelDecayMode(){return vv_pdfs_channelDecayMode;};
		const	vector< vector<int> > & Get_vv_productionMode(){return vv_productionMode;};
		const	vector< vector<int> > & Get_vv_pdfs_productionMode(){return vv_pdfs_productionMode;};
		int DecayMode(const std::string & s);
		int DecayModeFromProcessName(const std::string & s);
		int ProductionMode(const std::string & s);
		void Set_Cv_i(int i){_Cv_i=i;};
		void Set_Cf_i(int i){_Cf_i=i;};
		int Get_Cv_i(){return _Cv_i;};
		int Get_Cf_i(){return _Cf_i;};
		double ScaleCvCfHiggs(int countingOrParametric, int dm, int pm, int c, int s, double bs, const double *par);// counting 1,  parametric 2
		SMHiggsBuilder* GetSMHiggsBuilder(){if(_smhb)return _smhb; else {cout<<"ERROR SMHiggsBuilder not inited yet !"<<endl;exit(1);}};
		void SetFlatPars(double *pars);
		double CalcGammaTot();
		vector<structPOI> AddCvCf(vector<TString> scv, vector<TString> scf);
		
			const vector< vector<double> > & Get_v_Pars(){return v_Pars;} // on all flat parameters (nuisances and pois)  ,  index is the same as v_pdftype
		vector<structPOI> AddCX(vector<TString> scv, vector<TString> scg, vector<TString> sct, vector<TString> scb, vector<TString> scgl);// X  can be 4, 5, 7 ...
		double ScaleCXHiggs(int countingOrParametric, int dm, int pm, int c, int s, double bs, const double *par);// counting 1,  parametric 2
		int Get_Cgg_i(){return _Cgg_i;};
		int Get_Cvv_i(){return _Cvv_i;};
		int Get_Ctt_i(){return _Ctt_i;};
		int Get_Cbb_i(){return _Cbb_i;};
		int Get_Cglgl_i(){return _Cglgl_i;};

		int Get_par_i(TString spar){for(int i=0; i<v_uncname.size(); i++) if(v_uncname[i]==spar)return i+1;};

		TString GetErrEstAlgo(){return _ErrEstAlgo;};
		void SetErrEstAlgo(TString s){_ErrEstAlgo=s;};
		void ShowCvCfHiggsScales(double *par);
		void Set_printPdfEvlCycle(double d){_printPdfEvlCycle=d;};
		void PrintParametricChannelDataEntries();



		//  /**********  upgrade for TH1 based input   **********/
		public:
			void AddChannel(string channel_name, vector<string> vprocname , vector<TH1D*> sigTHs, 
					vector<TH1D*> bkgTHs, int decaymode=-1);
			// TH  norm unc:  LogNormal and TruncatedGaussian 
			void AddUncertaintyTH(string chname, int index_sample, double uncertainty_in_relative_fraction_down, double uncertainty_in_relative_fraction_up, int pdf_type, std::string uncname );
			// TH  norm unc:  From SideBand
			void AddUncertaintyTH(string chname, int index_sample, double rho, double rho_err, double B, int pdf_type, std::string uncname );
			// TH  shape parameters
			void AddUncertaintyTH(string chname, int index_sample, vector<TH1D*> par, int pdf_type, std::string uncname );
			// set data
			void AddObservedDataTH(int index_channel, TH1D* th);
			void AddObservedDataTH(string c, TH1D* th);
			void SetDataTH(VChannelTH  data, bool bRealData=true){v_data_th=data; if(bRealData)v_data_real_th=data;}
			void SetProcessNamesTH(int ch, vector<std::string> vproc); //need make check
			void SetProcessNamesTH(string ch, vector<std::string> vproc); //need make check
			TH1D* GetExpectedTH(string ch, string proc);
	
			
			const VChannelVSampleVUncertaintyVParameterTH& Get_vvvv_uncpar_th(){return vvvv_uncpar_th;}
			void Set_vvvv_uncpar_th(const VChannelVSampleVUncertaintyVParameterTH& vvvv){vvvv_uncpar_th=vvvv;}
			const VChannelVSampleVUncertainty& Get_vvv_pdftype_th(){return vvv_pdftype_th;}
			const VChannelVSampleVUncertainty& Get_vvv_idcorrl_th(){return vvv_idcorrl_th;}

			const VChannelVSampleTH& Get_vv_exp_sigbkgs_nonscaled_th() {return vv_exp_sigbkgs_th;}
			const VChannelVSampleTH& Get_vv_exp_sigbkgs_th()           {return vv_exp_sigbkgs_scaled_th;}
			int NumOfHistChannels(){return vv_exp_sigbkgs_th.size();}
			int GetNSigprocInChannelTH(int i){return v_sigproc_th.at(i);} //need make check
			const vector< bool > & Get_v_statusUpdated_th(){return v_statusUpdated_th;};
			const vector<int>& GetListOfShapeUncertaintiesTH(int c, int p){ return vvv_shapeuncindex_th[c][p]; } // need make check

			const vector<std::string>& GetProcessNamesTH(int ch){return vv_procname_th[ch];}  //need make check
			std::string GetChannelNameTH(int i){return v_channelname_th.at(i);} //need make check 
			const	vector< vector<int> > & Get_vv_channelDecayModeTH(){return vv_channelDecayMode_th;};
			const	vector< vector<int> > & Get_vv_productionModeTH(){return vv_productionMode_th;};
			const VChannelTH& Get_v_dataTH(){return v_data_th;}
		private:
			vector<int> v_sigproc_th;
			vector<std::string> v_channelname_th; // start from 0 
			vector< vector<std::string> > vv_procname_th; //name of each process in each channel
			VChannelTH v_data_th; // could be pseudo-data for bands
			VChannelTH v_data_real_th; // real data, not changed during entire run 
			VChannelVSampleTH vv_exp_sigbkgs_th;
			VChannelVSampleTH vv_exp_sigbkgs_scaled_th;
			VChannelVSampleTH vv_sigbkgs_varied_th;
			VChannelTH v_data_asimovb_th; // asimov data b only
			VChannelTH v_data_asimovsb_th; // asimov data mu*s + b 
			VChannelVSampleTH vv_randomized_sigbkgs_th;   
			VChannelVSampleTH vv_randomized_sigbkgs_scaled_th;
			VChannelVSampleTH vv_fitted_sigbkgs_th;    // fitted in data with mu=0
			VChannelVSampleTH vv_fitted_sigbkgs_scaled_th;   // fitted in data with mu being tested
			VChannelVSampleVUncertaintyVParameterTH  vvvv_uncpar_th;
			VChannelVSampleVUncertainty vvv_pdftype_th;
			VChannelVSampleVUncertainty vvv_idcorrl_th;
			vector< vector<int> > vv_channelDecayMode_th;
			vector< vector<int> > vv_productionMode_th;
			vector< vector< vector<int> > > vvv_shapeuncindex_th; //channel:process:shapeunc
			vector< vector<std::pair<int, int> > > vvp_th_connectNuisBinProc;// keep in memory:  a nuisance affects a list of [channel, process]
			vector< bool > v_statusUpdated_th;  //// monitoring if nuisances belonging to it(each pdf/process) updated 


		private:
			VDChannel v_data; // could be pseudo-data for bands
			VDChannel v_data_real; // real data, not changed during entire run 
			VDChannel v_data_asimovb; // asimov data b only
			VDChannel v_data_asimovsb; // asimov data mu*s + b 
			VChannelVSample vv_exp_sigbkgs;
			VChannelVSample vv_exp_sigbkgs_scaled;
			VChannelVSample vv_randomized_sigbkgs;   
			VChannelVSample vv_randomized_sigbkgs_scaled;
			VChannelVSample vv_fitted_sigbkgs;    // fitted in data with mu=0
			VChannelVSample vv_fitted_sigbkgs_scaled;   // fitted in data with mu being tested
			VChannelVSampleVUncertaintyVParameter  vvvv_uncpar;
			VChannelVSampleVUncertainty vvv_pdftype;
			VChannelVSampleVUncertainty vvv_idcorrl;

			vector<std::string> v_channelname; // start from 0 

			CRandom *_rdm;
			bool b_systematics;
			bool b_ForceSymmetryError;
			bool b_MultiSigProcShareSamePDF;

			vector<double> v_TruncatedGaussian_maxUnc;// record the maximum uncertainty for each uncorrelated source
			vector<int> v_pdftype; // 0th = -1,  if take a pdftype for idcorrl, then the indice is idcorrl
			double _common_signal_strength;
			int max_uncorrelation; // 
			void ConfigUncertaintyPdfs();

			bool b_AllowNegativeSignalStrength;

			vector<double> v_GammaN; // record the number of sideband(or MC raw) events for each uncorrelated source

			MapStrV map_flatPars;  // mainly on normalization
			vector< vector<double> > v_Pars; // on all flat parameters (nuisances and pois)  ,  index is the same as v_pdftype
			// for flat:   0. modified value,  1. init value,  2., 3., 4,  
			vector<int>  v_flatparId; // only store the flat nuisances which are not in workspaces 

			vector<std::string> v_uncname; // start from 0,   if take a name for idcorrl, then the indice is idcorrl-1; 
			vector<bool> v_uncFloatInFit; // start from 0,   if take a name for idcorrl, then the indice is idcorrl-1; 
			int _debug;

			vector<int> v_sigproc; // number of signal processes in each channel
			vector< vector<std::string> > vv_procname; //name of each process in each channel

			std::string _modelName;
			int _decayMode;

			vector< vector< vector<int> > > vvv_shapeuncindex; //channel:process:shapeunc

			bool bMoveUpShapeUncertainties;

			double * _norminalPars; 
			double * _norminalParsSigma; 
			double * _randomizedPars; 

			double * _fittedParsInData_bonly; // perform a fit on real data  with mu=0	
			double * _fittedParsInData_sb;	  // perform a fit on real data  with mu being tested
			double * _fittedParsInData_global;	  // perform a fit on real data  with mu floating
			double * _fittedParsInPseudoData_bonly;	// perform a fit on pseudo data with mu=0
			double * _fittedParsInPseudoData_sb;	// perform a fit on pseudo data with mu being tested
			double * _fittedParsInPseudoData_global;	// perform a fit on pseudo data with mu floating
			int _tossToyConvention;
		       	// convention 0 for the LEP type we used to do;
		        // convention 1 for the LHC type agreed on LHC-HCG meeting https://indico.cern.ch/getFile.py/access?contribId=48&sessionId=5&resId=0&materialId=slides&confId=139132
			int _UseBestEstimateToCalcQ; // this should sit inside CLsBase....   however to avoid a global CLsBase, we put it inside CountingModel   FIXME  low priority 
			// 0 use randomized measurements;  1 use the original best estimates;  2 use the fitted sets in data

			//// start to add unbinned parametric shape into model,  for only 1-dimention 
			bool bRedefineObservableRange; int ObservableBins; double ObservableRangeMin, ObservableRangeMax;
			bool bHasParametricShape;// for both binned and unbinned 
			vector< vector<string> > vv_pdfs; // every process has its own pdf, in each parametric channel
			vector< vector<double> > vv_pdfs_norm; //normalization of each pdf, i.e. expected number of events in each process 
			vector< vector<double> > vv_pdfs_norm_scaled; //normalization of each pdf, i.e. expected number of events in each process 
			vector< vector<double> > vv_pdfs_norm_randomized; //normalization of each pdf, i.e. expected number of events in each process 
			vector< vector<double> > vv_pdfs_norm_randomized_scaled; //normalization of each pdf, i.e. expected number of events in each process 
			vector< vector<double> > vv_pdfs_norm_fitted; //normalization of each pdf, i.e. expected number of events in each process    fitted in data with mu=0
			vector< vector<double> > vv_pdfs_norm_fitted_scaled; //normalization of each pdf, i.e. expected number of events in each process   fitted in data with mu being tested
			// uncertainties ....   each source affects parameters  -->  two additional sets of parameters
			vector< vector< vector<int> > > vvv_pdfs_idcorrl;
			vector< vector< vector<int> > > vvv_pdfs_pdftype;
			vector< vector< vector< vector<double> > > > vvv_pdfs_normvariation; // correponding to normalization changes of the uncertainty source
			vector< vector<double> > vv_pdfs_norm_varied; //normalization of each pdf, i.e. expected number of events in each process 
			vector<int> v_pdfs_sigproc;// number of signal processes in each channel
			vector<string> v_pdfs_channelname;

			vector< TString > v_pdfs_sb; // tot pdf for s+b model    in each channel 
			vector< TString > v_pdfs_b; //  tot pdf for b-only model  in each channel
			vector< TString > v_pdfs_s; //  tot pdf for s-only model  in each channel

			vector<TString> v_pdfs_observables; // observable in each channel
			vector<RooAbsData*>  v_pdfs_roodataset_toy;
			vector<RooAbsData*>  v_pdfs_roodataset; // could be pseudo-data for bands
			vector<RooAbsData*>  v_pdfs_roodataset_real; // real data, not changed during entire run 
			vector<RooAbsData*>  v_pdfs_roodataset_tmp;
			vector<RooAbsData*>  v_pdfs_roodataset_asimovb; // b-only  asimov dataset
			vector<RooAbsData*>  v_pdfs_roodataset_asimovsb; // b-only  asimov dataset


			vector< vector<TString> > vv_pdfs_normNAME;
			vector< vector<TString> > vv_pdfs_extranormNAME;

			RooWorkspace * _w;
			RooWorkspace * _w_varied;

			vector<int> v_pdfs_floatParamsIndcorr;      // only for params
			vector<int> v_pdfs_floatParamsType;      // only for params
			vector<string> v_pdfs_floatParamsName;     // only for params
			vector<double> v_pdfs_floatParamsVaried;  // only for params
			vector< vector<double> > v_pdfs_floatParamsUnc; // from 0 to max_uncorl


			// following are not supported yet 
			vector< vector<int> > vv_pdfs_npar; // number of parameters for each pdf
			vector< vector<double *> > vv_pdfs_params; //nominal parameters for each pdf
			vector< vector< vector<double*> > > vvv_pdfs_params_up; // correponding to up shift of the uncertainty source, e.g. jet energy scale
			vector< vector< vector<double*> > > vvv_pdfs_params_down; //  correponding to down shift of the uncertainty source, e.g. jet energy scale
			vector< vector<double *> > vv_pdfs_params_varied; //nominal parameters for each pdf
			vector< vector< vector<double*> > > vvv_pdfs_params_min; // correponding to max shift of down side of the uncertainty source, e.g. jet energy scale
			vector< vector< vector<double*> > > vvv_pdfs_params_max; //  correponding to max shift of up side of the uncertainty source, e.g. jet energy scale
			vector< vector< vector<RooRealVar*> > > vvv_pdfs_paramsRRV; //nominal parameters for each pdf

			vector< int > v_pdfs_nbin; // need for throwing random numbers,  should consistent for all pdfs in a channel, by default 100 
			vector< double > v_pdfs_xmin;
			vector< double > v_pdfs_xmax;
			vector< vector<double> > vv_pdfs_data; // in each channel, it has a list of events 
			vector< vector<double> > vv_pdfs_data_toy; // in each channel, it has a list of events
			// three types of uncertainties: 1. only affect shape;  2. only affect normalization; 3. affect both 
			vector< vector< vector<int> > > vvv_pdfs_unctype;
			vector< vector< vector<int> > > vvv_pdfs_nuisancesindex;

			vector< vector< vector< vector<double> > > > vvvv_pdfs_ChProcSetEvtVals;
			vector< vector< vector< vector<double> > > > vvvv_pdfs_ChProcSetParVals;
			vector< vector<int> > vv_pdfs_curSetIndex;

			vector< vector<string> > vv_pdfs_procname;

			MapStrVV map_param_sources; // map <paramName, vector< idcorrl, sigmaL, sigmaR > > 


			vector< vector< bool > > vv_pdfs_statusUpdated;// monitoring if nuisances belonging to it(each pdf/process) updated
			vector< bool > v_pdfs_statusUpdated;// monitoring if nuisances belonging to it(each pdf/process) updated
			vector< vector< bool > > vv_statusUpdated;  // 
			vector< vector<std::pair<int, int> > > vvp_pdfs_connectNuisBinProc;// keep in memory:  a nuisance affects a list of [channel, process]
			vector< vector<std::pair<int, int> > > vvp_connectNuisBinProc;// keep in memory:  a nuisance affects a list of [channel, process]

			vector< vector< int > > TMP_vvpdfs_chprocINT; 


			bool _bFixingPOIs;
			vector< std::pair<TString, double> > vsPOIsToBeFixed;
			
			vector<structPOI> vPOIs; 
			

			int minuitSTRATEGY;
			double minuitTolerance;
			double nuisancesRange;
			int maximumFunctionCallsInAFit;
			int maxsets_forcaching;
			int PrintParameterFrom;
			int PrintParameterTo;

			int minuitPrintLevel;

			int _PhysicsModel;

			bool noErrorEstimation; // running minos, but without error estimation


			vector< vector<int> > vv_channelDecayMode;
			vector< vector<int> > vv_pdfs_channelDecayMode;
			vector< vector<int> >vv_productionMode;
			vector< vector<int> >vv_pdfs_productionMode;

			int _Cv_i, _Cf_i;  // for CvCfHiggs

			SMHiggsBuilder *_smhb;
			double _HiggsMass; // when no "MH" parameter, then set it to model 
			double _GammaTot; // for recalc total decay width when coupling changes 
			double _CvCf_gg; // for recalc total decay width when coupling changes 
			double _CvCf_zg; // for recalc total decay width when coupling changes 
			double *_pardm; // to replace vrdm in FluctuatedNumbers

			int _MH_i;

			int _Cvv_i, _Cbb_i, _Ctt_i, _Cgg_i, _Cglgl_i, _Ctptp_i, _Czgzg_i, _Cff_i;

			TString _ErrEstAlgo;
			double _printPdfEvlCycle;

	};
	CountingModel* CombineModels(CountingModel *cms1, CountingModel *cms2);
}
#endif   /* ----- #ifndef CountingModel_H----- */
