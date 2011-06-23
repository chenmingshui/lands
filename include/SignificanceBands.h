#ifndef  SIGNIFICANCEBANDS_H
#define  SIGNIFICANCEBANDS_H
#include "CLsLimit.h"
namespace lands{
	class SignificanceBands{
		public:
			SignificanceBands(CLsBase *cb, CountingModel* cms);
			~SignificanceBands();
			void Bands(int noutcome);
			double GetSignificance(int i){return _rcls[i+2];} // -2, -1, 0, 1, 2
			double GetSignificanceMean(){return _rcls[5];}

			void SetDebug(int d) {_debug=d;}
			void SetModel(CountingModel* cms){_cms=cms;}

			const vector<double>& GetDifferentialSignificances(){return _difrcls;}
		private:
			void Bands();
			int _noutcomes;
			vector<double> _difrcls, _vrcls, _vpcls; 
			int _debug;
			CLsBase *_frequentist;
			CountingModel *_cms;
			double _rcls[6];
			
	};
}
#endif   /* ----- #ifndef SIGNIFICANCEBANDS_H  ----- */

