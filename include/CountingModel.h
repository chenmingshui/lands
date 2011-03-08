#ifndef CountingModel_H
#define CountingModel_H
/*
 * =====================================================================================
 *
 *       Filename: CountingModel.h 
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  03/10/2010 10:50:40 PM CET
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Mingshui Chen (), Mingshui.Chen@cern.ch
 *        Company:  Univsity of Florida
 *
 * =====================================================================================
 */
#include <vector>
#include <string>
#include "CRandom.h"
#include "PdfRandom.h"
#include "Utilities.h"

namespace lands{
	typedef vector<double> VDChannel;
	//typedef vector<int> VIChannel;
	typedef vector<double> VIChannel;
	typedef vector< vector<double> > VChannelVSample;
	typedef vector< vector< vector<int> > > VChannelVSampleVUncertainty;
	typedef vector< vector< vector< vector<double> > > > VChannelVSampleVUncertaintyVParameter;

	enum enumPdfType {typeLogNormal=1, typeTruncatedGaussian=2, typeControlSampleInferredLogNormal=11 };
	
	class CountingModel
	{ 
		public:
			CountingModel();
			~CountingModel();
			void AddChannel(std::string channel_name, double num_expected_signal, double num_expected_bkg_1, double num_expected_bkg_2=-1, 
						double num_expected_bkg_3=-1, double num_expected_bkg_4=-1, double num_expected_bkg_5=-1, double num_expected_bkg_6 = -1 );
			void AddChannel(double num_expected_signal, double num_expected_bkg_1, double num_expected_bkg_2=-1, 
						double num_expected_bkg_3=-1, double num_expected_bkg_4=-1, double num_expected_bkg_5=-1, double num_expected_bkg_6 = -1 );

			void AddChannel(double num_expected_signal, vector<double> num_expected_bkgs);
			void AddChannel(std::string channel_name, double num_expected_signal, vector<double> num_expected_bkgs);
			void AddChannel(vector<double> num_expected_yields);
			void AddChannel(std::string channel_name, vector<double> num_expected_yields);

			// LogNormal and TruncatedGaussian 
			void AddUncertainty(int index_channel, int index_sample, double uncertainty_in_relative_fraction, int pdf_type, int index_correlation );

			// From SideBand
			void AddUncertainty(int index_channel, int index_sample, double rho, double rho_err, double B, int pdf_type, int index_correlation );
			// when B is large, it's close to Gaussian. then use Gaussian, don't use Gamma function


			void AddObservedData(int index_channel, double num_data);
			void SetData(VDChannel data){v_data=data;};
			void SetData(vector<int> data){for(int i=0; i<data.size(); i++) v_data[i]=data[i];};

			VChannelVSample Get_vv_exp_sigbkgs(){return vv_exp_sigbkgs_scaled;};
			VDChannel Get_v_data(){return v_data;};

			VDChannel Get_v_exp_sigbkgs(int channel){return vv_exp_sigbkgs_scaled[channel];};

			double GetExpectedNumber(int index_channel, int index_sample);

			void Print(int printLevel=0);
			bool Check();

			VChannelVSample FluctuatedNumbers();	
			VIChannel GetToyData_H0();// background only hypothesis
			VIChannel GetToyData_H1();// alternative hypothesis
			
			int NumOfChannels(){return vv_exp_sigbkgs.size();};
			void SetUseSystematicErrors(bool b){b_systematics=b;ConfigUncertaintyPdfs();};
			bool IsUsingSystematicsErrors(){return b_systematics;};

			void SetSignalScaleFactor(double r);
			double GetSignalScaleFactor(){return _common_signal_strength; };
			void SetRdm(CRandom *rdm){_rdm=rdm;};
			CRandom * GetRdm(){return _rdm;};

			void UseAsimovData(int b=0);

			CountingModel CombineModels(CountingModel *cms1, CountingModel *cms2);

			VChannelVSampleVUncertaintyVParameter Get_vvvv_uncpar(){return vvvv_uncpar;};
			VChannelVSampleVUncertainty Get_vvv_pdftype(){return vvv_pdftype;};
			VChannelVSampleVUncertainty Get_vvv_idcorrl(){return vvv_idcorrl;};

			void RemoveChannelsWithExpectedSignal0orBkg0(int king = 2); // king=0 for removing a bin with bkg=0, 1 for sig=0, 2 for bkg=0||sig=0  

			int Get_max_uncorrelation() {return max_uncorrelation;};
			VDChannel Get_v_TruncatedGaussian_maxUnc() {return v_TruncatedGaussian_maxUnc;};
			vector<int> Get_v_pdftype() {return v_pdftype;};
		private:
			VDChannel v_data;
			VChannelVSample vv_exp_sigbkgs;
			VChannelVSample vv_exp_sigbkgs_scaled;
			VChannelVSampleVUncertaintyVParameter  vvvv_uncpar;
			VChannelVSampleVUncertainty vvv_pdftype;
			VChannelVSampleVUncertainty vvv_idcorrl;

			vector<std::string> v_channelname;

			CRandom *_rdm;
			bool b_systematics;

			vector<PdfRandom*> v_TruncatedGaussian;
			vector<double> v_TruncatedGaussian_maxUnc;
			vector<int> v_pdftype;
			double _common_signal_strength;
			int max_uncorrelation;
			void ConfigUncertaintyPdfs();

	};
};
#endif   /* ----- #ifndef CountingModel_H----- */
