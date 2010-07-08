/*
 *    Description:  doing frequentist approach + CLs(not fully frequentist)
 *			w/ and w/o systematic/statistic errors 
 */

#ifndef  CLSLIMIT_H
#define  CLSLIMIT_H
#include "CRandom.h"
#include <vector>
#include "CountingModel.h"
using std::vector;
using std::pair;
namespace lands{
	class CLsBase
	{
		public:
			CLsBase();
			~CLsBase();
			void SetModel(CountingModel *model){_model=model;};
			bool BuildM2lnQ(int ntoys=100000, int sbANDb_bOnly_sbOnly=0);
			bool BuildM2lnQ(CountingModel *model, int ntoys=100000, int sbANDb_bOnly_sbOnly=0);// 0 for sbANDb, 1 for bOnly, 2 for sbOnly
			void SetRdm(CRandom *rdm);
			vector<double> Get_m2logQ_b();
			vector<double> Get_m2logQ_sb();
			double Get_m2lnQ_data();

			// data, observed 
			double CLsb();double CLs();double CLb();
			double CLb(double lnq );
			double PValue(double lnq );
			//expected  bkg only 
			double CLsb_b();double CLs_b();	double CLb_b();

			void tmpFun0(vector<double> & vlogQ, vector<double>& vlogQ_prob);
			double lnQsb_sigma(int sigma); //-2, -1, 0, 1, 2, 3 // 0-median, 3-mean
			void CheckFractionAtHighEnd(vector<double> vlogQ, vector<double> vlogQ_prob);
			double SignificanceComputation(int ntoys_for_sb, int ntoys_for_b);
			double SignificanceComputation(int ntoys_for_sb, int ntoys_for_b, vector<double>& vsignificance, vector<double> & vsignificance_cp);
			double SignificanceForData(int ntoys_for_b);


			void SetDebug(int debug);

			CRandom* GetRdm();

			void SetLogQ_b(vector<double> vlnQ_b);	
			void SetLogQ_sb(vector<double> vlnQ_sb);	
			void SetLogQ_data(double lnQ_data);	

		private:
			void ProcessM2lnQ();
			double *Q_b; double *Q_sb;
			int *iq_b; // sort m2logQ_b
			int *iq_sb; // sort m2logQ_sb
			double Q_b_exp; double Q_b_data;// double Q_b_median; double Q_b_mean;
			double _nsig;	double _nbkg; 	double _ndat; 
			int _nexps;
			int _nchannels;
			CRandom *_rdm;
			int _debug;
			CountingModel *_model;

	};
	class CLsLimit
	{
		public:
			CLsLimit(){_debug=0; _alpha = 0.05; _clstolerance=0.001; };
			~CLsLimit(){}; 	
			void SetAlpha(double alpha); // Confidence Level = 1 - alpha

			double LimitOnSignalScaleFactor(CountingModel *cms, CLsBase *frequentist, int nexps=100000);
			double LimitOnSignalScaleFactor(CountingModel *cms, double minRtoScan, double maxRtoScan, CLsBase *frequentist, int nexps=100000, int nsteps = 10);
			double GetLimit();
			vector<double> GetvTestedScaleFactors(); //
			vector<double> GetvTestedCLs();//	

			void DoingStatisticalBandsForLimit(CountingModel *cms, CLsBase *frequentist, int nexps_to_buildM2lnQ=100000, int nexps_to_doStatisticalBands=1000);
			double Limit_sigma(int nsigma); // return limit at -2sigma, -1sigma, median(50%), 1sigma, 2sigma
			double Limit_mean(); //return average value mathmatically.....
			vector<double> GetDifferentialLimits();
			vector<double> GetvLimits(); // corresponding to all possible outcomes  ,  cummulative
			vector<double> GetvLimits_CP(); // corresponding to all possible outcomes

			void DoingStatisticalBandsForCLs(CountingModel *cms, CLsBase *frequentist, int nexps_to_buildM2lnQ = 100000);
			double CLs_sigma(int nsigma); // return CLs at -2sigma, -1sigma, median(50%), 1sigma, 2sigma
			double CLs_mean(); //return average value mathmatically.....
			vector<double> GetDifferentialCLsReq1();
			vector<double> GetvCLsReq1(); // corresponding to all possible outcomes
			vector<double> GetvCLsReq1_CP(); // corresponding to all possible outcomes

			void SetDebug(int debug);
			CLsBase* GetFrequentist();

			void SetCLsTolerance(double tolerance = 0.001 );


		private:
			vector<double> _vR;
			vector<double> _vCLs; // for a fixed s,b,d, trying to converge at CLs=0.05 to get r=r95%, filling all CLs produced during that process
			double _r95;

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
	};

};
#endif   /* ----- #ifndef CLsLIMIT_H  ----- */

