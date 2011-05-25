#ifndef  LIMITBANDS_H
#define  LIMITBANDS_H
#include "CLsLimit.h"
#include "BayesianBase.h"
namespace lands{
	class LimitBands{
		public:
			LimitBands();
			LimitBands(CLsLimit *cl, CLsBase *cb, CountingModel* cms);
			LimitBands(BayesianBase *bb, CountingModel* cms);
			LimitBands(CLsLimit *cl, CLsBase *cb, BayesianBase *bb, CountingModel* cms);
			~LimitBands();
			void Bands(double alpha, int noutcome, bool doCLs, int ntoysM2lnQ, bool doBys, int ntoysBys);
			void CLsLimitBands(double alpha=0.05, int noutcome=1000, int ntoysM2lnQ=100000);
			void BysLimitBands(double alpha=0.05, int noutcome=1000, int ntoysBys=20000);
			double GetCLsLimit(int i){return _rcls[i+2];}; // -2, -1, 0, 1, 2
			double GetBysLimit(int i){return _rbys[i+2];}; // -2, -1, 0, 1, 2
			double GetCLsLimitMean(){return _rcls[5];};
			double GetBysLimitMean(){return _rbys[5];};

			void SetDebug(int d) {_debug=d;};
			void SetAlpha(double a){fAlpha=a; fConfidenceLevel-fAlpha;};
			void SetCLs(CLsLimit* clslimit, CLsBase* frequentist) {_clslimit=clslimit; _frequentist=frequentist;};
			void SetBys(BayesianBase * bys ){_byslimit=bys;};
			void SetModel(CountingModel* cms){_cms=cms;};
			void SetTossPseudoDataConvention(int i){_tossPseudoDataConvention=i;};

			vector<double> GetDifferentialLimitsCLs(){return _difrcls;};
			vector<double> GetDifferentialLimitsBys(){return _difrbys;};
		private:
			void Bands();
			bool _doCLs, _doBys;
			double fAlpha;
			double fConfidenceLevel;
			int _noutcomes;
			int _actualOutComesForBys;
			int _ntoysM2lnQ;
			int _ntoysBys;
			vector<double> _difrcls, _vrcls, _vpcls; 
			vector<double> _difrbys, _vrbys, _vpbys; 
			int _debug;
			CLsBase *_frequentist;
			CLsLimit *_clslimit;
			BayesianBase *_byslimit;
			CountingModel *_cms;
			double _rcls[6], _rbys[6];
			int _tossPseudoDataConvention;
			
	};
};
#endif   /* ----- #ifndef LIMITBANDS_H  ----- */

