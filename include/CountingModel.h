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
#include "RooAbsPdf.h"
#include "RooRealVar.h"
#include "RooBifurGauss.h"
#include "RooWorkspace.h"
#include <map>
#include "TString.h"

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

	enum enumPdfType {typeLogNormal=1, typeTruncatedGaussian=2, typeGamma=3, typeShapeGaussianLinearMorph=4, typeShapeGaussianQuadraticMorph=5, 
		typeBifurcatedGaussian=6, typeFlat=7, typeControlSampleInferredLogNormal=11 };
	// for uncertainties only affecting Shape, we can choose different morphing algorithm.  in commom, user must provide three templates:  norminal,  shift_1sigma_up, shift_1sigma_down
	// i.e. 3 parameters for each shape uncertainty in each bin .... ,  the interface alway do normalization (to unity) for all three templates.
	
	enum enumPhysicsModel {typeStandardModel=1, typeChargedHiggs=2};

	struct structPOI {
		TString name;
		double value;
		double errUp;
		double errDown;
		structPOI(TString s, double v, double eu, double ed):name(s),value(v),errUp(eu),errDown(ed){}
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
			void AddChannel(std::string channel_name, vector<double> num_expected_yields, int signal_processes=1);
			void AddChannel(std::string channel_name, vector<double> num_expected_signals, vector<double> num_expected_bkgs);

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

			void SetProcessNames(int ch, vector<std::string> vproc); //need make check
			void SetProcessNames(string ch, vector<std::string> vproc); //need make check
			const vector<std::string>& GetProcessNames(int ch){return vv_procname[ch];}  //need make check
			void SetModelName(const std::string& s){_modelName = s;}
			const std::string& GetModelName(){return _modelName;}
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
			void ConstructAsimovData(int b, bool nominal=true);
			const vector< RooAbsData* >& Get_v_pdfs_roodataset_asimovsb(){return v_pdfs_roodataset_asimovsb;} // in each channel, it has a list of events
			const vector< double >& Get_v_pdfs_floatParamsVaried(){return v_pdfs_floatParamsVaried;}
			const vector< vector< double > >& Get_v_pdfs_floatParamsUnc(){ return v_pdfs_floatParamsUnc;}
			void Set_v_pdfs_floatParamsUnc(const vector< vector< double > >& vv ){ v_pdfs_floatParamsUnc = vv;}
			void SetFlatParameterRange(int id, double middle, double errLow, double errUp);
			void SetFlatNormalizationRange(int id, double errLow, double errUp);
			const vector<int>& Get_v_pdfs_floatParamsIndcorr() {return v_pdfs_floatParamsIndcorr;}      // only for params
			const vector<int>& Get_v_pdfs_floatParamsType() {return v_pdfs_floatParamsType;}      // only for params

			const vector<string>& Get_v_pdfs_channelname(){return v_pdfs_channelname;}
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

			void AddChannel(string channel_name, RooRealVar* observable, vector<RooAbsPdf*> sigPdfs, vector<double> sigNorms, vector<RooAbsArg*> vsExtraNorm, 
				       	vector<RooAbsPdf*> bkgPdfs, vector<double> bkgNorms, vector<RooAbsArg*> vbExtraNorm);
			// need to add names of each parameter .... 
			double EvaluateLnQ(int ch, int dataOrToy); // for Likelihood ratio
			double EvaluateChi2(double *par, vector< vector< vector<float> > >& vvv_cachPdfValues, int bUseBestEstimateToCalcQ=1);          // for Chi2
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

			void SetPhysicsModel(int i){_PhysicsModel=i; if(i!=1 && i!=2){ cout<<"ERROR: we only support typeStandardModel and typeChargedHiggs"<<endl; exit(1);}};
			int GetPhysicsModel(){return _PhysicsModel;};

			void ForceSymmetryError(bool b){b_ForceSymmetryError = b;};
			void MultiSigProcShareSamePDF(bool b){b_MultiSigProcShareSamePDF = b;};


			void DumpFitResults(double *pars, TString ssave);
			void Set_maxSetsCaching(int i){maxsets_forcaching=i;};
			void Set_PrintParameter(int i, int j){PrintParameterFrom =i; PrintParameterTo = j;};
			int GetPrintParameterFrom(){return PrintParameterFrom;};
			int GetPrintParameterTo(){return PrintParameterTo;};


			bool bFixingPOIs(){return _bFixingPOIs;};
			const vector< std::pair<TString, double> > & GetPOIsToBeFixed(){return vsPOIsToBeFixed;}
			void SetbFixingPOIs(bool b){_bFixingPOIs = b; if(!b) vsPOIsToBeFixed.clear();};
			void SetPOItoBeFixed(TString sp, double p){ vsPOIsToBeFixed.push_back(std::make_pair(sp,p)); _bFixingPOIs = true;};
			const vector<structPOI> & POIs(){return vPOIs;};
			void addPOI(structPOI p) {vPOIs.push_back(p);};
			void setPOI(int i, double v, double eu, double ed){vPOIs[i].value=v; vPOIs[i].errUp=eu; vPOIs[i].errDown=ed;};

			void Set_minuitPrintLevel(int i){minuitPrintLevel = i;};
			int Get_minuitPrintLevel(){return minuitPrintLevel;};
		private:
			VDChannel v_data; // could be pseudo-data for bands
			VDChannel v_data_real; // real data, not changed during entire run 
			VDChannel v_data_asimovb; // asimov data b only
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

			vector<std::string> v_uncname; // start from 0,   if take a name for idcorrl, then the indice is idcorrl-1; 
			vector<bool> v_uncFloatInFit; // start from 0,   if take a name for idcorrl, then the indice is idcorrl-1; 
			int _debug;

			vector<int> v_sigproc; // number of signal processes in each channel
			vector< vector<std::string> > vv_procname; //name of each process in each channel

			std::string _modelName;

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

	};
	CountingModel* CombineModels(CountingModel *cms1, CountingModel *cms2);
}
#endif   /* ----- #ifndef CountingModel_H----- */
