#include "BinnedInterface.h"
#include "TH1D.h"
#include <iostream>
using namespace std;
namespace lands{
	void BinnedInterface::AddChannel(TH1* h_signal, TH1* h_background1, TH1* h_background2, TH1* h_background3, TH1* h_background4, TH1* h_background5){
		
		_hsignals.push_back(h_signal);
		// nbins of all samples should be identical, and min/max of xaxis shoud also be identical for all smaples, need check here, FIXME
		int nbins = h_signal->GetXaxis()->GetNbins();	
//		cout<<"nbins="<<nbins<<endl;
		for(int i=1; i<=nbins; i++){  // do we want to account underflow/overflow bins ?
//			cout<<"bin"<<i<<"= "<<h_signal->GetBinContent(i)<<endl;
			_cms->AddChannel(
				h_signal->GetBinContent(i),
				h_background1->GetBinContent(i),
				h_background2?(h_background2->GetBinContent(i)):-1,
				h_background3?(h_background3->GetBinContent(i)):-1,
				h_background4?(h_background4->GetBinContent(i)):-1,
				h_background5?(h_background5->GetBinContent(i)):-1
					);
		}
	}	
	void BinnedInterface::AddData(int index_channel, TH1* h_data){
		// need check if h_data bins are consistent with signal and background // FIXME
		int nbins=_hsignals.at(index_channel)->GetXaxis()->GetNbins();	
		int startchannel = 0;
		for(int ich = 0; ich<index_channel; ich++){
			startchannel+=_hsignals.at(ich)->GetXaxis()->GetNbins();		
		}
		for( int i = 1 ; i <=nbins; i++){
//			cout<<"AddData bin "<<i<<"="<<h_data->GetBinContent(i)<<endl;
			_cms->AddObservedData(i+startchannel-1, h_data->GetBinContent(i) );
		}
	}
	void BinnedInterface::AddUncertainty(int index_channel, int index_sample, double error_in_relative_fraction, int pdf_type, int index_correlation){
		// need check if h_data bins are consistent with signal and background // FIXME
		int nbins=_hsignals.at(index_channel)->GetXaxis()->GetNbins();	
		int startchannel = 0;
		for(int ich = 0; ich<index_channel; ich++){
			startchannel+=_hsignals.at(ich)->GetXaxis()->GetNbins();		
		}
		for( int i = 1 ; i <=nbins; i++){
//			cout<<"AddUncertainty bin "<<i<<"="<<error_in_relative_fraction<<endl;
			_cms->AddUncertainty(i+startchannel-1, index_sample, error_in_relative_fraction, pdf_type, index_correlation );
		}

	}
	void BinnedInterface::AddUncertainty(int index_channel, int index_sample, TH1* h_error_in_relative_fraction, int pdf_type, int index_correlation){
		// need check if h_data bins are consistent with signal and background // FIXME
		int nbins=_hsignals.at(index_channel)->GetXaxis()->GetNbins();	
		int startchannel = 0;
		for(int ich = 0; ich<index_channel; ich++){
			startchannel+=_hsignals.at(ich)->GetXaxis()->GetNbins();		
		}
		for( int i = 1 ; i <=nbins; i++)
			_cms->AddUncertainty(i+startchannel-1, index_sample, h_error_in_relative_fraction->GetBinContent(i), pdf_type, index_correlation );
	}
};

