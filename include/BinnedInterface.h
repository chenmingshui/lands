#ifndef  BINNEDINTERFACE_H
#define  BINNEDINTERFACE_H
#include "TH1.h"
#include "CountingModel.h"
#include <vector>
namespace lands{
	class BinnedInterface{
		public:
			BinnedInterface(CountingModel *cms){_cms=cms; _hsignals.clear(); };
			~BinnedInterface(){_cms=0;_hsignals.clear();};
			void AddChannel(TH1* h_signal, TH1* h_background1, TH1* h_background2=0, TH1* h_background3=0, TH1* h_background4=0, TH1* h_background5=0);	
			void AddData(int index_channel, TH1* h_data);
			void AddUncertainty(int index_channel, int index_sample, double error_in_relative_fraction, int pdf_type, int index_correlation);
			void AddUncertainty(int index_channel, int index_sample, TH1* h_error_in_relative_fraction, int pdf_type, int index_correlation);
			CountingModel* GetModel(){return _cms;};	
		private:
			CountingModel *_cms;
			vector<TH1*> _hsignals;
	};
};
#endif   /* ----- #ifndef BINNEDINTERFACE_INC  ----- */

