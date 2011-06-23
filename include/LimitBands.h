#ifndef  LIMITBANDS_H
#define  LIMITBANDS_H
#include "CLsLimit.h"
#include "BayesianBase.h"
#include "TString.h"
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
			double GetCLsLimit(int i){return _rcls[i+2];} // -2, -1, 0, 1, 2
			double GetBysLimit(int i){return _rbys[i+2];} // -2, -1, 0, 1, 2
			double GetCLsLimitMean(){return _rcls[5];}
			double GetBysLimitMean(){return _rbys[5];}

			void SetDebug(int d) {_debug=d;}
			void SetAlpha(double a){fAlpha=a; fConfidenceLevel-fAlpha;}
			void SetCLs(CLsLimit* clslimit, CLsBase* frequentist) {_clslimit=clslimit; _frequentist=frequentist;}
			void SetBys(BayesianBase * bys ){_byslimit=bys;}
			void SetModel(CountingModel* cms){_cms=cms;}
			void SetTossPseudoDataConvention(int i){_tossPseudoDataConvention=i;}

			const vector<double>& GetDifferentialLimitsCLs(){return _difrcls;}
			const vector<double>& GetDifferentialLimitsBys(){return _difrbys;}
			void IsM2lnQGridPreComputed(bool b, TString s){bM2lnQGridPreComputed=b; sFileM2lnQGrid=s;}
			void SetPlotLevel(int i){_plotLevel = i;}

			void Set_bOnlyEvalCL_forVR(bool b){bOnlyEvalCL_forVR=b;}
			void Set_vR_toEval(const vector<double>& v){_vR_toEval = v;}
			const vector<vector<double> >& Get_vvCL_forVR(){return _vvCL_forVR;}
			const vector<double>& Get_vR_toEval(){return _vR_toEval;}

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
			bool bM2lnQGridPreComputed;
			TString sFileM2lnQGrid;

			int _plotLevel;

			vector<double> _vR_toEval;
			bool bOnlyEvalCL_forVR;
			vector< vector<double> > _vvCL_forVR;
			
	};
}
#endif   /* ----- #ifndef LIMITBANDS_H  ----- */

