#ifndef BAYESIANBASE_H
#define BAYESIANBASE_H
#include <vector>
#include "CountingModel.h"
using namespace std;
namespace lands{
	typedef enum { //prior on the cross section 
		flat = 10, // uncorrelated  // default
		corr = 20, // correlated  
		prior_1overSqrtS = 30 // 1/sqrt(s)
	} PRIOR;
	class BayesianBase{
		public:
			BayesianBase();
			BayesianBase(double alpha, double precision);
			BayesianBase(CountingModel *cms, double alpha=0.05, double precision=1.e-3);
			~BayesianBase();

			double Limit(CountingModel *cms);
			double Limit(double alpha=0.05);
			void SetModel(CountingModel *cms);

			//single channel counting experiment, limit in average which is with mathmatical meaning
			double LimitMeanInMath(); 

			void SetAlpha(double alpha){fAlpha=alpha;fConfidenceLevel=1-fAlpha;};
			void SetDebug(int debug){_debug=debug;};
			void SetNumToys(int ntoys) {_nexps_to_averageout_sys=ntoys;};
			void GenToys();
			double AverageIntegral(double rlow);
			double glintegral(double xlow, int iexps); 
			double GetAlpha(){return fAlpha;};
			double GetLimit(){return _limit;};
			double Likelihood(double r);
			double PosteriorPdf(int bins=100, double rmin=0, double rmax=0); 
			// will return computational error on the r due to iteration cutout, the fPrecision,  in case of no systematics
			double ErrorOnR_DueToPrecision();
			double ErrorOnR_DueToFiniteToys(int ntrials, double& mean);
			vector<double> GetVR(){return _vr;};
			vector<double> GetVP(){return _vp;};
			void SetCrossSectionPrior(PRIOR prior){_prior=prior;};			
			PRIOR GetCrossSectionPrior(){return _prior;};			
		private:
			double fPrecision;
			double fAlpha;
			double fConfidenceLevel; // = 1- fAlpha;
			int _debug;
			CountingModel *_cms;
			int _nexps_to_averageout_sys;

			int _nchannels;

			vector<double *> _vs, _vb; // new double[_nexps_to_averageout_sys][_nchannels];
			double *_d; // new double[_nchannels];

			// Gauss-Laguerre Quadrature   // from "genlimit" package, CDF note 7587
			int _ngl;
			double *_xgl, *_lwgl;
			double _logscale;
			double _norm;
			double *_stot, *_btot;
			double _limit;
			vector<double> _vr, _vp;
			double _dtot;
			PRIOR _prior;
	};
};

#endif //#ifndef BAYESIANBASE_H
