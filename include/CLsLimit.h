#ifndef  CLSLIMIT_H
#define  CLSLIMIT_H
#include "CRandom.h"
#include <vector>
#include "CountingModel.h"
#include "TMinuit.h"
#include "TGraphErrors.h"
/*
 *    Description:  doing frequentist approach + CLs(not fully frequentist)
 *			w/ and w/o systematic/statistic errors 
 */

using std::vector;
using std::pair;
namespace lands{
	extern CountingModel * cms_global;
	extern TMinuit * myMinuit;
	extern vector<double>  vdata_global;
	extern TF1* fitToRvsCL_expo;
	extern double _significances;
	extern double *_inputNuisances;
	extern double *_inputNuisancesSigma;
	extern double *_startNuisances; // _xxxxNuisances[0] is the POI signal strength
	extern double *_minNuisances;
	extern double *_maxNuisances;
	extern double _customRMin;
	extern double _customRMax;
	extern bool _bPositiveSignalStrength;
	extern vector<double> _lastParams;
	extern vector<double> _currParams;
	extern vector< vector< vector<float> > >  vvv_cachPdfValues; // cached pdf values for layers Channel<Process<Event>>>
	extern vector< vector< double > >  vv_cachCountingParts; // cached likelihood values for layers Channel<Process>
	extern double _countPdfEvaluation;
	extern bool _bDumpFinalFitResults;
	void Chisquare(Int_t &npar, Double_t *gin, Double_t &f,  Double_t *par, Int_t iflag); 
	double MinuitFit(int model, double &r, double &er, double mu=0, double *pars=0, bool hasBestFitted = false, int debug=0, int *success=0);
	bool DoAfit(double mu, vector<double> vdata, vector<RooAbsData*> vrds, double* pars);
	class CLsBase
	{
		public:
			CLsBase();
			~CLsBase();
			void SetModel(CountingModel *model){_model=model; cms_global = model;}
			bool BuildM2lnQ(int ntoys=100000, int sbANDb_bOnly_sbOnly=0, bool reUsePreviousToys=false);
			bool BuildM2lnQ_b(int ntoys=100000, bool reUsePreviousToys=false, bool bWriteToys = false);
			bool BuildM2lnQ_sb(int ntoys=100000, bool reUsePreviousToys=false, bool bWriteToys = false);
			bool BuildM2lnQ_data();
			bool BuildM2lnQ(CountingModel *model, int ntoys=100000, int sbANDb_bOnly_sbOnly=0, bool reUsePreviousToys=false);// 0 for sbANDb, 1 for bOnly, 2 for sbOnly
			void SetRdm(CRandom *rdm);
			vector<double> Get_m2logQ_b();
			vector<double> Get_m2logQ_sb();
			double Get_m2lnQ_data();

			// data, observed 
			double CLsb(double &err);double CLs(double & err);double CLb(double &err);
			double CLb(double lnq, double &err );
			double PValue(double lnq );
			//expected  bkg only 
			double CLsb_b();double CLs_b();	double CLb_b();

			void tmpFun0(vector<double> & vlogQ, vector<double>& vlogQ_prob);
			double lnQsb_sigma(int sigma); //-2, -1, 0, 1, 2, 3 // 0-median, 3-mean
			void CheckFractionAtHighEnd(vector<double> vlogQ, vector<double> vlogQ_prob);
			double SignificanceComputation(int ntoys_for_sb, int ntoys_for_b);
			double SignificanceComputation(int ntoys_for_sb, int ntoys_for_b, vector<double>& vsignificance, vector<double> & vsignificance_cp);
			double SignificanceForData(int ntoys_for_b);
			double SignificanceForData(double qdata, vector<double> vq);


			void SetDebug(int debug);

			CRandom* GetRdm();

			void SetLogQ_b(vector<double> vlnQ_b);	
			void SetLogQ_sb(vector<double> vlnQ_sb);	
			void SetLogQ_data(double lnQ_data);	

			void SetTestStatistics(int ts = 1);
			int GetTestStatistics(){return test_statistics;}

			int GetNexps(){return _nexps;}

			void printM2LnQInfo(int sbANDb_bOnly_sbOnly);
			void checkFittedParsInData(bool bReadPars=false, bool bWritePars=false, TString sfilename="");
			double M2lnQ(bool& success, int checkFailure=0, int dataOrToy=1);
			void prepareLogNoverB();

			double FindLimitFromTGE(TGraphErrors *tge, double alpha, double &limit, double & limitErr, TString plotName="");
			double FindLimitFromPreComputedGrid(std::map<double, TTree*> gridCLsb, std::map<double, TTree*> gridCLb, std::map<double, double> gridQdata, double alpha, TString plotName=""); // from precomputed m2lnQ grid to extract r corresponding to _alpha ... e.g. 0.05


		private:
			void ProcessM2lnQ();
			double *Q_b; double *Q_sb;
			int *iq_b; // sort m2logQ_b
			int *iq_sb; // sort m2logQ_sb
			double Q_b_exp; double Q_b_data;// double Q_b_median; double Q_b_mean;
			double _nsig;	double _nbkg; 	double _ndat; 
			double *_lognoverb;
			int _nexps;
			int _nchannels;
			CRandom *_rdm;
			int _debug;
			CountingModel *_model;

			int test_statistics;  // 1 for Q_LEP,  2 for Q_TEV, 3 for Q_ATLAS,  5 for the LHC type

	};
	class CLsLimit
	{
		public:
			CLsLimit(){_debug=0; _alpha = 0.05; _clstolerance=0.001; _rule = 1; bAdaptiveSampling=false; fAdditionalNToysFactor=1.;}  // by default, we use CLs instead of CLsb
			~CLsLimit(){} 	
			void SetAlpha(double alpha); // Confidence Level = 1 - alpha

			double LimitOnSignalScaleFactor(CountingModel *cms, CLsBase *frequentist, int nexps=100000);
			double LimitOnSignalScaleFactor(CountingModel *cms, double minRtoScan, double maxRtoScan, CLsBase *frequentist, int nexps=100000, int nsteps = 10);
			double GetLimit();
			const vector<double>& GetvTestedScaleFactors(); //
			const vector<double>& GetvTestedCLs();//	

			void DoingStatisticalBandsForLimit(CountingModel *cms, CLsBase *frequentist, int nexps_to_buildM2lnQ=100000, int nexps_to_doStatisticalBands=1000);
			double Limit_sigma(int nsigma); // return limit at -2sigma, -1sigma, median(50%), 1sigma, 2sigma
			double Limit_mean(); //return average value mathmatically.....
			const vector<double>& GetDifferentialLimits();
			const vector<double>& GetvLimits(); // corresponding to all possible outcomes  ,  cummulative
			const vector<double>& GetvLimits_CP(); // corresponding to all possible outcomes

			void DoingStatisticalBandsForCLs(vector<double> vsb, vector<double> vb);
			double CLs_sigma(int nsigma); // return CLs at -2sigma, -1sigma, median(50%), 1sigma, 2sigma
			double CLs_mean(); //return average value mathmatically.....
			const vector<double>& GetDifferentialCLsReq1();
			const vector<double>& GetvCLsReq1(); // corresponding to all possible outcomes
			const vector<double>& GetvCLsReq1_CP(); // corresponding to all possible outcomes

			void SetDebug(int debug);
			CLsBase* GetFrequentist();

			void SetCLsTolerance(double tolerance = 0.001 );
			void SetRule(int rule = 1);

			double LimitErr(){return _r95err;}


			double FeldmanCousins(CountingModel *cms, double minRtoScan, double maxRtoScan, CLsBase *frequentist, int nexps=100000, int nsteps = 10);
			const vector< vector<double> >& GetFCconstruction(){return _FCconstruction;}
			void SetAdaptiveSampling(bool b){bAdaptiveSampling=b;}
			void SetAdditionalNToysFactor(double d){fAdditionalNToysFactor=d;}
			double * Get_CLsProjected(){return _CLsProjected;};

		private:
			vector<double> _vR;
			vector<double> _vCLs; // for a fixed s,b,d, trying to converge at CLs=0.05 to get r=r95%, filling all CLs produced during that process
			double _r95;
			double _r95err;

			int _nexps;  // for everytime to cal CLs
			int _npossibleoutcomes;// do projecting
			CLsBase *_frequentist;

			vector<double> _vR95, _vR95_CP;  // cummulative
			vector<double> _vCLs_Req1, _vCLs_Req1_CP; //cummulative
			vector<double> _differentialR95s; 
			vector<double> _differentialCLs_Req1; 
			double _r95Projected[6]; // -2, -1, 0(median), 1, 2, mean 
			double _CLsProjected[6]; // -2, -1, 0(median), 1, 2, mean 
			int _debug;
			double _alpha;// Confidence Level = 1 - alpha
			double _clstolerance;

			int _rule;  // 1 for CLs, 2 for CLsb

			// in each vector<float>, first increasing order of -2lnQ,   last element is the r being tested, followed by  q_up and q_data
			vector< vector<double> > _FCconstruction;

			// implemented for FeldmanCousin limit  setting
			bool bAdaptiveSampling; // doing sampling with adaptive number of toys
			double fAdditionalNToysFactor; // allow user to require more toys when doing adaptive sampling. 
	};

}
#endif   /* ----- #ifndef CLsLIMIT_H  ----- */

