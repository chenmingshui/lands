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
			void SetSignalStrength(int toysAtDifferentSignalStrength, double injectingSignalStrength){
				_toysAtDifferentSignalStrength = toysAtDifferentSignalStrength; 
				_injectingSignalStrength = injectingSignalStrength;	
			}
		private:
			void Bands();
			int _noutcomes;
			vector<double> _vrcls, _vpcls; //_difrcls, 
			int _debug;
			//CLsBase *_frequentist;
			CountingModel *_cms;
			double _rcls[6];
			double _injectingSignalStrength ; 
			int _toysAtDifferentSignalStrength ; // b-only,    1 to use injectingSignalStrength,  2 to use best fitted signal strength	
	};
}
#endif   /* ----- #ifndef MLLXSBANDS_H  ----- */

