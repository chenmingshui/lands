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

using namespace std;
namespace lands{
	typedef vector<double> VDChannel;
	//typedef vector<int> VIChannel;
	typedef vector<double> VIChannel;
	typedef vector< vector<double> > VChannelVSample;
	typedef vector< vector< vector<int> > > VChannelVSampleVUncertainty;
	typedef vector< vector< vector< vector<double> > > > VChannelVSampleVUncertaintyVParameter;

	enum enumPdfType {typeLogNormal=1, typeTruncatedGaussian=2, typeGamma=3, typeShapeGaussianLinearMorph=4, typeShapeGaussianQuadraticMorph=5, 
		typeBifurcatedGaussian=6, typeControlSampleInferredLogNormal=11 };
	// for uncertainties only affecting Shape, we can choose different morphing algorithm.  in commom, user must provide three templates:  norminal,  shift_1sigma_up, shift_1sigma_down
	// i.e. 3 parameters for each shape uncertainty in each bin .... ,  the interface alway do normalization (to unity) for all three templates.
	
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
			void SetData(VDChannel data, bool bRealData=true){v_data=data; if(bRealData)v_data_real=data;};
			void SetData(TString filename, TString datahistname, bool bRealData=true);
			void SetData(vector<int> data){for(int i=0; i<data.size(); i++) v_data[i]=data[i];};

			VChannelVSample Get_vv_exp_sigbkgs_nonscaled(){return vv_exp_sigbkgs;};
			VChannelVSample Get_vv_exp_sigbkgs(){return vv_exp_sigbkgs_scaled;};
			void Set_vv_randomized_sigbkgs(VChannelVSample vv){vv_randomized_sigbkgs=vv;};
			VChannelVSample Get_vv_randomized_sigbkgs(){return vv_randomized_sigbkgs_scaled;};
			void Set_vv_fitted_sigbkgs(VChannelVSample vv){vv_fitted_sigbkgs=vv;};
			VChannelVSample Get_vv_fitted_sigbkgs(){return vv_fitted_sigbkgs;};
			void Set_vv_fitted_sigbkgs_scaled(VChannelVSample vv){vv_fitted_sigbkgs_scaled=vv;};
			VChannelVSample Get_vv_fitted_sigbkgs_scaled(){return vv_fitted_sigbkgs_scaled;};
			VDChannel Get_v_data(){return v_data;};
			VDChannel Get_v_data_real(){return v_data_real;};

			VDChannel Get_v_exp_sigbkgs(int channel){return vv_exp_sigbkgs_scaled[channel];};

			double GetExpectedNumber(int index_channel, int index_sample);

			void Print(int printLevel=0);
			bool Check();

			VChannelVSample FluctuatedNumbers(double *pars = NULL, bool scaled=true, int bUseBestEstimateToCalcQ=1); 
			// 0 use randomized set,  1 use best estimates a priori, 2 use fitted posterior (after looking at data)
			
			VIChannel GetToyData_H0(double *pars=NULL);// background only hypothesis
			VIChannel GetToyData_H1(double *pars=NULL);// alternative hypothesis
			
			int NumOfChannels(){return vv_exp_sigbkgs.size();};
			void SetUseSystematicErrors(bool b){b_systematics=b;ConfigUncertaintyPdfs();};
			bool IsUsingSystematicsErrors(){return b_systematics;};

			void SetSignalScaleFactor(double r, bool bScaleBestEstimate=true);
			double GetSignalScaleFactor(){return _common_signal_strength; };
			void SetRdm(CRandom *rdm){_rdm=rdm;};
			CRandom * GetRdm(){return _rdm;};

			void UseAsimovData(int b=0);


			VChannelVSampleVUncertaintyVParameter Get_vvvv_uncpar(){return vvvv_uncpar;};
			VChannelVSampleVUncertainty Get_vvv_pdftype(){return vvv_pdftype;};
			VChannelVSampleVUncertainty Get_vvv_idcorrl(){return vvv_idcorrl;};

			int RemoveChannelsWithExpectedSignal0orBkg0(int king = 2); // king=0 for removing a bin with bkg=0, 1 for sig=0, 2 for bkg=0||sig=0  

			int Get_max_uncorrelation() {return max_uncorrelation;};
			VDChannel Get_v_TruncatedGaussian_maxUnc() {return v_TruncatedGaussian_maxUnc;};
			vector<int> Get_v_pdftype() {return v_pdftype;};

			void SetAllowNegativeSignalStrength(bool b){b_AllowNegativeSignalStrength = b;};
			bool AllowNegativeSignalStrength(){return b_AllowNegativeSignalStrength;};

			vector<double> Get_v_GammaN(){return v_GammaN;};
			vector<std::string> Get_v_uncname(){return v_uncname;};
			void SetDebug(int i){_debug=i;};
			int GetDebug(){return _debug;};
			std::string GetChannelName(int i){return v_channelname.at(i);}; //need make check 
			int GetNSigprocInChannel(int i){return v_sigproc.at(i);}; //need make check
			vector<int> Get_v_sigproc(){return v_sigproc;};
			void SetProcessNames(int ch, vector<std::string> vproc); //need make check
			void SetProcessNames(string ch, vector<std::string> vproc); //need make check
			vector<std::string> GetProcessNames(int ch){return vv_procname[ch];};  //need make check
			void SetModelName(std::string s){_modelName = s;};
			std::string GetModelName(){return _modelName;};
			void MakeListOfShapeUncertainties();
			vector<int> GetListOfShapeUncertainties(int c, int p){ return vvv_shapeuncindex[c][p]; }; // need make check
			void SetMoveUpShapeUncertainties(bool b){bMoveUpShapeUncertainties=b;};
			bool GetMoveUpShapeUncertainties(){return bMoveUpShapeUncertainties;};

			void SetTossToyConvention(int c){_tossToyConvention = c;};
			int GetTossToyConvention(){return _tossToyConvention;};
			double* Get_fittedParsInData_sb(){return _fittedParsInData_sb;};
			double* Get_fittedParsInData_b(){return _fittedParsInData_bonly;};
			void Set_fittedParsInData_sb(double *p){
			       	_fittedParsInData_sb=p;
				vv_fitted_sigbkgs_scaled = FluctuatedNumbers(p);
				vv_pdfs_norm_fitted_scaled = vv_pdfs_norm_varied;
			};
			void Set_fittedParsInData_b(double *p){
				_fittedParsInData_bonly=p;
				vv_fitted_sigbkgs= FluctuatedNumbers(p, 0);
				vv_pdfs_norm_fitted = vv_pdfs_norm_varied;
			};

			void SetUseBestEstimateToCalcQ(bool b=true){ _UseBestEstimateToCalcQ=b;};
			bool UseBestEstimateToCalcQ(){ return _UseBestEstimateToCalcQ;};

			// start to add parametric shape into model
			bool hasParametricShape(){return bHasParametricShape;};
			vector< vector<string> > Get_vv_pdfs(){return vv_pdfs;};
			vector< vector<double *> > Get_vv_pdfs_params(){return vv_pdfs_params;}; //nominal parameters for each pdf
			vector< vector<double> > Get_vv_pdfs_norm_nonscaled(){return vv_pdfs_norm;}; //normalization of each pdf, i.e. expected number of events in each process 
			vector< vector<int> > Get_vv_pdfs_npar(){return vv_pdfs_npar;}; // number of parameters for each pdf
			vector< int > Get_v_pdfs_nbin(){return v_pdfs_nbin;}; // need for throwing random numbers,  should consistent for all pdfs in a channel, by default 100 
			vector< double > Get_v_pdfs_xmin(){return v_pdfs_xmin;};
			vector< double > Get_v_pdfs_xmax(){return v_pdfs_xmax;};
			vector< vector<double> > Get_vv_pdfs_data(){return vv_pdfs_data;}; // in each channel, it has a list of events
			// uncertainties ....   each source affects parameters  -->  two additional sets of parameters
			vector< vector< vector<int> > > Get_vvv_pdfs_idcorrl(){return vvv_pdfs_idcorrl;};
			vector< vector< vector<int> > > Get_vvv_pdfs_pdftype(){return vvv_pdfs_pdftype;};
			// three types of uncertainties: 1. only affect shape;  2. only affect normalization; 3. affect both 
			vector< vector< vector<int> > > Get_vvv_pdfs_unctype(){return vvv_pdfs_unctype;};
			vector< vector< vector<double*> > > Get_vvv_pdfs_params_up(){return vvv_pdfs_params_up;}; // correponding to up shift of the uncertainty source, e.g. jet energy scale
			vector< vector< vector<double*> > > Get_vvv_pdfs_params_down(){return vvv_pdfs_params_down;}; //  correponding to down shift of the uncertainty source, e.g. jet energy scale
			vector< vector< vector< vector<double> > > > Get_vvv_pdfs_normvariation(){return vvv_pdfs_normvariation;}; // correponding to normalization changes of the uncertainty source
			vector< vector<double *> > Get_vv_pdfs_params_varied(){return vv_pdfs_params_varied;}; //nominal parameters for each pdf
			vector< vector<double> > Get_vv_pdfs_norm_varied(){return vv_pdfs_norm_varied;}; //normalization of each pdf, i.e. expected number of events in each process 
			vector< vector<double> > Get_vv_pdfs_norm_scaled(){return vv_pdfs_norm_scaled;}; //normalization of each pdf, i.e. expected number of events in each process 
			vector< vector<double> > Get_vv_pdfs_data_toy(){return vv_pdfs_data_toy;}; // in each channel, it has a list of events
			void Set_vv_pdfs_norm_randomized(VChannelVSample vv){vv_pdfs_norm_randomized=vv;};
			VChannelVSample Get_vv_pdfs_norm_randomized(){return vv_pdfs_norm_randomized_scaled;};
			void Set_vv_pdfs_norm_fitted(VChannelVSample vv){vv_pdfs_norm_fitted=vv;};
			VChannelVSample Get_vv_pdfs_norm_fitted(){return vv_pdfs_norm_fitted;};
			void Set_vv_pdfs_norm_fitted_scaled(VChannelVSample vv){vv_pdfs_norm_fitted_scaled=vv;};
			VChannelVSample Get_vv_pdfs_norm_fitted_scaled(){return vv_pdfs_norm_fitted_scaled;};
			vector<int> Get_v_pdfs_sigproc(){return v_pdfs_sigproc;};
			vector< RooDataSet* > Get_v_pdfs_roodataset_toy(){return v_pdfs_roodataset_toy;}; // in each channel, it has a list of events
			vector< RooDataSet* > Get_v_pdfs_roodataset(){return v_pdfs_roodataset;}; // in each channel, it has a list of events
			vector< RooDataSet* > Get_v_pdfs_roodataset_real(){return v_pdfs_roodataset_real;}; // in each channel, it has a list of events
			vector< double > Get_v_pdfs_floatParamsVaried(){return v_pdfs_floatParamsVaried;};
			vector< vector< double > >Get_v_pdfs_floatParamsUnc(){ return v_pdfs_floatParamsUnc;};
			vector<int> Get_v_pdfs_floatParamsIndcorr() {return v_pdfs_floatParamsIndcorr;};      // only for params

			vector<string> Get_v_pdfs_channelname(){return v_pdfs_channelname;};
			vector<TString> Get_v_pdfs_sb(){return v_pdfs_sb;};
			vector<TString> Get_v_pdfs_s(){return v_pdfs_s;};
			vector<TString> Get_v_pdfs_b(){return v_pdfs_b;};
			vector<TString> Get_v_pdfs_observables(){return v_pdfs_observables;};
			RooWorkspace * GetWorkSpace(){return _w;};


			void AddChannel(string channel_name, RooRealVar* observable, vector<RooAbsPdf*> sigPdfs, vector<double> sigNorms, vector<RooAbsPdf*> bkgPdfs, vector<double> bkgNorms, RooWorkspace *w );
			// need to add names of each parameter .... 
			double EvaluateLnQ(int ch, int dataOrToy); // for Likelihood ratio
			double EvaluateChi2(double *par, bool bUseBestEstimateToCalcQ=true);          // for Chi2
			double EvaluateGL(int ch, double xr); // for bayesian 
			double EvaluateGL(vector< vector<double> > vvnorms, vector<double> vparams, double xr); // for bayesian 
			void AddObservedDataSet(int index_channel, RooDataSet* rds);
			void AddObservedDataSet(string channelname, RooDataSet* rds);
			void SetDataForUnbinned(vector< RooDataSet*> data, bool bRealData=true);
			void SetDataForUnbinned(TString filename, bool bRealData=true);
			void SetTmpDataForUnbinned(vector< RooDataSet*> data); // set it before doing fit
			void SetToyForUnbinned(vector< RooDataSet*> data); // set it before doing fit
			void AddUncertaintyOnShapeNorm(int index_channel, int index_sample, double uncertainty_in_relative_fraction_down, double uncertainty_in_relative_fraction_up, int pdf_type, int index_correlation );
			void AddUncertaintyOnShapeNorm(int index_channel, int index_sample, double uncertainty_in_relative_fraction_down, double uncertainty_in_relative_fraction_up, int pdf_type, string uncname);
			void AddUncertaintyOnShapeNorm(string chname, int index_sample, double uncertainty_in_relative_fraction_down, double uncertainty_in_relative_fraction_up, int pdf_type, string uncname);
			void AddUncertaintyOnShapeParam(string pname, double mean, double sigmaL, double sigmaR, double rangeMin=0, double rangeMax=0 );

			void AddUncertaintyAffectingShapeParam(string uname, string pname, double mean, double sigmaL, double sigmaR, double rangeMin, double rangeMax );
		private:
			VDChannel v_data; // could be pseudo-data for bands
			VDChannel v_data_real; // real data, not changed during entire run 
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

			vector<double> v_TruncatedGaussian_maxUnc;// record the maximum uncertainty for each uncorrelated source
			vector<int> v_pdftype; // 0th = -1,  if take a pdftype for idcorrl, then the indice is idcorrl
			double _common_signal_strength;
			int max_uncorrelation; // 
			void ConfigUncertaintyPdfs();

			bool b_AllowNegativeSignalStrength;

			vector<double> v_GammaN; // record the number of sideband(or MC raw) events for each uncorrelated source

			vector<std::string> v_uncname; // start from 0,   if take a name for idcorrl, then the indice is idcorrl-1; 
			int _debug;

			vector<int> v_sigproc; // number of signal processes in each channel
			vector< vector<std::string> > vv_procname; //name of each process in each channel

			std::string _modelName;

			vector< vector< vector<int> > > vvv_shapeuncindex; //channel:process:shapeunc

			bool bMoveUpShapeUncertainties;

			double * _fittedParsInData_bonly; // perform a fit on real data  with mu=0	
			double * _fittedParsInData_sb;	  // perform a fit on real data  with mu being tested
			double * _fittedParsInData_global;	  // perform a fit on real data  with mu floating
			double * _fittedParsInPseudoData_bonly;	// perform a fit on pseudo data with mu=0
			double * _fittedParsInPseudoData_sb;	// perform a fit on pseudo data with mu being tested
			double * _fittedParsInPseudoData_global;	// perform a fit on pseudo data with mu floating
			int _tossToyConvention;
		       	// convention 0 for the LEP type we used to do;
		        // convention 1 for the LHC type agreed on LHC-HCG meeting https://indico.cern.ch/getFile.py/access?contribId=48&sessionId=5&resId=0&materialId=slides&confId=139132
			bool _UseBestEstimateToCalcQ; // this should sit inside CLsBase....   however to avoid a global CLsBase, we put it inside CountingModel   FIXME  low priority 

			//// start to add unbinned parametric shape into model,  for only 1-dimention 
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
			vector<RooDataSet*>  v_pdfs_roodataset_toy;
			vector<RooDataSet*>  v_pdfs_roodataset; // could be pseudo-data for bands
			vector<RooDataSet*>  v_pdfs_roodataset_real; // real data, not changed during entire run 
			vector<RooDataSet*>  v_pdfs_roodataset_tmp;


			vector< vector<TString> > vv_pdfs_normNAME;

			RooWorkspace * _w;
			RooWorkspace * _w_varied;

			vector<int> v_pdfs_floatParamsIndcorr;      // only for params
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
			vector< vector<string> > vv_pdfs_procname;
	};
	CountingModel* CombineModels(CountingModel *cms1, CountingModel *cms2);
};
#endif   /* ----- #ifndef CountingModel_H----- */
