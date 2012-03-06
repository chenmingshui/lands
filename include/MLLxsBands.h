#ifndef  MLLXSBANDS_H
#define  MLLXSBANDS_H
#include "CLsLimit.h"
namespace lands{
	class MLLxsBands{
		public:
			MLLxsBands(CountingModel* cms);
			~MLLxsBands();
			void Bands(int noutcome);
			double GetMLLxs(int i){return _rcls[i+2];} // -2, -1, 0, 1, 2
			double GetMLLxsMean(){return _rcls[5];}

			void SetDebug(int d) {_debug=d;}
			void SetModel(CountingModel* cms){_cms=cms;}

			const vector<double>& GetDifferentialMLLxs(){return _vrcls;}
		private:
			void Bands();
			int _noutcomes;
			vector<double> _vrcls, _vpcls; //_difrcls, 
			int _debug;
			//CLsBase *_frequentist;
			CountingModel *_cms;
			double _rcls[6];
			
	};
}
#endif   /* ----- #ifndef MLLXSBANDS_H  ----- */

