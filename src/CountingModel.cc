/*
 * =====================================================================================
 *
 *       Filename: CountingModel.cc 
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  03/10/2010 10:52:47 PM CET
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Mingshui Chen (), Mingshui.Chen@cern.ch
 *        Company:  Univsity of Florida
 *
 * =====================================================================================
 */
#include "CountingModel.h"
#include <iostream>
#include <cstdio>
#include <cmath>
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "RooAbsPdf.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooBifurGauss.h"
#include "RooWorkspace.h"
#include "RooStats/RooStatsUtils.h"

using namespace std;
using namespace RooFit;
using namespace RooStats;
namespace lands{
	CountingModel::CountingModel(){
		v_data.clear();
		vv_exp_sigbkgs.clear();
		vv_exp_sigbkgs_scaled.clear();
		vvvv_uncpar.clear();
		vvv_pdftype.clear();
		vvv_idcorrl.clear();
		v_channelname.clear();
		_rdm=0;
		b_systematics=0;
		v_TruncatedGaussian_maxUnc.clear();
		v_pdftype.clear();
		_common_signal_strength=1;
		max_uncorrelation=0;
		b_AllowNegativeSignalStrength = 1;
		v_GammaN.clear();
		v_uncname.clear();
		v_sigproc.clear();
		vv_procname.clear();
		_debug=0;
		_modelName = "model";
		vvv_shapeuncindex.clear();
		bMoveUpShapeUncertainties = false;

		bHasParametricShape = false;
		vv_pdfs.clear();	
		vv_pdfs_params.clear();	
		vv_pdfs_norm.clear();	
		vv_pdfs_norm_scaled.clear();	
		vv_pdfs_npar.clear();	
		v_pdfs_nbin.clear();	
		v_pdfs_xmin.clear();	
		v_pdfs_xmax.clear();	
		vv_pdfs_data.clear();	
		vvv_pdfs_idcorrl.clear();
		vvv_pdfs_pdftype.clear();
		vvv_pdfs_unctype.clear();
		vvv_pdfs_params_up.clear();
		vvv_pdfs_params_down.clear();
		vvv_pdfs_normvariation.clear();
		vv_pdfs_params_varied.clear();
		vv_pdfs_norm_varied.clear();
		vv_pdfs_data_toy.clear();
		v_pdfs_sigproc.clear();
		v_pdfs_channelname.clear();
		vv_pdfs_procname.clear();
		v_pdfs_sb.clear();
		v_pdfs_b.clear();
		v_pdfs_s.clear();
		vvv_pdfs_params_min.clear();
		vvv_pdfs_params_max.clear();
		v_pdfs_observables.clear();
		v_pdfs_roodataset_toy.clear();
		v_pdfs_roodataset.clear();
		vv_pdfs_normNAME.clear();
		vvv_pdfs_paramsRRV.clear();

		v_pdfs_floatParamsName.clear();
		v_pdfs_floatParamsIndcorr.clear();
		v_pdfs_floatParamsUnc.clear();


		_w = new RooWorkspace();
		_w_varied = new RooWorkspace();
	}
	CountingModel::~CountingModel(){
		v_data.clear();
		vv_exp_sigbkgs.clear();
		vv_exp_sigbkgs_scaled.clear();
		vvvv_uncpar.clear();
		vvv_pdftype.clear();
		vvv_idcorrl.clear();
		v_channelname.clear();
		_rdm=0;

		v_TruncatedGaussian_maxUnc.clear();
		v_pdftype.clear();
		v_GammaN.clear();
		v_uncname.clear();
		v_sigproc.clear();
		vv_procname.clear();
		vvv_shapeuncindex.clear();

		//pdfs shapes
		for(int i=0; i<vv_pdfs.size(); i++){
		//	if(v_pdfs_roodataset[i]) delete v_pdfs_roodataset[i];
		//	if(v_pdfs_roodataset_toy[i]) delete v_pdfs_roodataset_toy[i];
		}
		vv_pdfs.clear();	
		vv_pdfs_params.clear();	
		vv_pdfs_norm.clear();	
		vv_pdfs_norm_scaled.clear();	
		vv_pdfs_npar.clear();	
		v_pdfs_nbin.clear();	
		v_pdfs_xmin.clear();	
		v_pdfs_xmax.clear();	
		vv_pdfs_data.clear();	
		vvv_pdfs_idcorrl.clear();
		vvv_pdfs_pdftype.clear();
		vvv_pdfs_unctype.clear();
		vvv_pdfs_params_up.clear();
		vvv_pdfs_params_down.clear();
		vvv_pdfs_normvariation.clear();
		vv_pdfs_params_varied.clear();
		vv_pdfs_norm_varied.clear();
		vv_pdfs_data_toy.clear();
		v_pdfs_sigproc.clear();
		v_pdfs_channelname.clear();
		vv_pdfs_procname.clear();
		v_pdfs_sb.clear();
		v_pdfs_b.clear();
		v_pdfs_s.clear();
		vvv_pdfs_params_min.clear();
		vvv_pdfs_params_max.clear();
		v_pdfs_observables.clear();
		v_pdfs_roodataset_toy.clear();
		v_pdfs_roodataset.clear();
		vv_pdfs_normNAME.clear();
		vvv_pdfs_paramsRRV.clear();
		v_pdfs_floatParamsName.clear();
		v_pdfs_floatParamsIndcorr.clear();
		v_pdfs_floatParamsUnc.clear();


	//	delete _w;
	//	delete _w_varied;
	}
	void CountingModel::AddChannel(std::string channel_name, double num_expected_signal, double num_expected_bkg_1, double num_expected_bkg_2, 
			double num_expected_bkg_3, double num_expected_bkg_4, double num_expected_bkg_5, double num_expected_bkg_6 ){
// FIXME need to work with bin content = 0; we add the channel anyway, then put it to be decided when calc  Q
		if( num_expected_signal <=0 ){
//			cout<<"You are adding a channel in which nsignal <= 0,   we will skip this channel"<<endl;
//			exit(0);
		}
		// should bkg_1 < 0, then nothing to do , you should put your bkgs as realistic
		if(num_expected_bkg_1 < 0 ) {cout<<"First Background < 0, exit "<<endl; exit(0);}
		if( 
				( num_expected_bkg_2 <0 && (num_expected_bkg_3 >=0 || num_expected_bkg_4>=0 || num_expected_bkg_5>=0 || num_expected_bkg_6>=0) ) ||
				( num_expected_bkg_3<0 && (num_expected_bkg_4>=0 ||num_expected_bkg_5>=0 || num_expected_bkg_6>=0) ) ||
				( num_expected_bkg_4 < 0 && ( num_expected_bkg_5>=0 || num_expected_bkg_6>=0 )  ) ||
				( num_expected_bkg_5 < 0 && num_expected_bkg_6>=0  )
		  ){
			cout<<"Something wrong with bkg inputs, exit "<<endl;
			 exit(0);
		}

		if( num_expected_bkg_1 <=0 &&
				(num_expected_bkg_2<=0 && num_expected_bkg_3 <=0 && num_expected_bkg_4<=0 && num_expected_bkg_5<=0 && num_expected_bkg_6<=0) ) {
		//	cout<<"Your are adding a channel with totbkg <= 0,  we will skip this channel"<<endl;
		//	exit(0);
		}

		if(channel_name==""){
			char tmp[256];
			sprintf(tmp, "channel_%d", v_channelname.size());
			channel_name==tmp;
			v_channelname.push_back(channel_name);

		}
		else v_channelname.push_back(channel_name);


		double tmp_totbkg = 0;
		
		vector<double> vsigbkgs; vsigbkgs.clear();
		vsigbkgs.push_back(num_expected_signal);
		vector< vector< vector<double> > > vvvuncpar; vvvuncpar.clear();
		vector< vector<int> > vvpdftype; vvpdftype.clear();
		vector< vector<int> > vvidcorrl; vvidcorrl.clear();

		vector< vector<double> > vvunc; vvunc.clear();
		vvvuncpar.push_back(vvunc);
		vector<int> vpdftype; vpdftype.clear();
		vvpdftype.push_back(vpdftype);
		vector<int> vidcorrl; vidcorrl.clear();
		vvidcorrl.push_back(vidcorrl);
		if(num_expected_bkg_1 >= 0 ){
			vsigbkgs.push_back(num_expected_bkg_1);
			vvvuncpar.push_back(vvunc);
			vvpdftype.push_back(vpdftype);
			vvidcorrl.push_back(vidcorrl);
			tmp_totbkg+=num_expected_bkg_1;
		}
		if(num_expected_bkg_2 >= 0 ){
			vsigbkgs.push_back(num_expected_bkg_2);
			vvvuncpar.push_back(vvunc);
			vvpdftype.push_back(vpdftype);
			vvidcorrl.push_back(vidcorrl);
			tmp_totbkg+=num_expected_bkg_2;
		}
		if(num_expected_bkg_3 >= 0 ){
			vsigbkgs.push_back(num_expected_bkg_3);
			vvvuncpar.push_back(vvunc);
			vvpdftype.push_back(vpdftype);
			vvidcorrl.push_back(vidcorrl);
			tmp_totbkg+=num_expected_bkg_3;
		}
		if(num_expected_bkg_4 >= 0 ){
			vsigbkgs.push_back(num_expected_bkg_4);
			vvvuncpar.push_back(vvunc);
			vvpdftype.push_back(vpdftype);
			vvidcorrl.push_back(vidcorrl);
			tmp_totbkg+=num_expected_bkg_4;
		}
		if(num_expected_bkg_5 >= 0 ){
			vsigbkgs.push_back(num_expected_bkg_5);
			vvvuncpar.push_back(vvunc);
			vvpdftype.push_back(vpdftype);
			vvidcorrl.push_back(vidcorrl);
			tmp_totbkg+=num_expected_bkg_5;
		}
		if(num_expected_bkg_6 >= 0 ){
			vsigbkgs.push_back(num_expected_bkg_6);
			vvvuncpar.push_back(vvunc);
			vvpdftype.push_back(vpdftype);
			vvidcorrl.push_back(vidcorrl);
			tmp_totbkg+=num_expected_bkg_6;
		}

		vv_exp_sigbkgs.push_back(vsigbkgs);
		vv_exp_sigbkgs_scaled.push_back(vsigbkgs);
		v_data.push_back(tmp_totbkg);	
		vvvv_uncpar.push_back(vvvuncpar);
		vvv_pdftype.push_back(vvpdftype);
		vvv_idcorrl.push_back(vvidcorrl);
		v_sigproc.push_back(1);
		vector<string> vproc; vproc.clear();
	       	for(int i=0; i<7; i++) {
			char tmp[256];
			sprintf(tmp, "%d", i);
			vproc.push_back(tmp);
		}
		vv_procname.push_back(vproc);
	}
	void CountingModel::AddChannel(double num_expected_signal, double num_expected_bkg_1, double num_expected_bkg_2, 
			double num_expected_bkg_3, double num_expected_bkg_4, double num_expected_bkg_5, double num_expected_bkg_6 ){
		AddChannel("",  num_expected_signal,  num_expected_bkg_1,  num_expected_bkg_2, 
				num_expected_bkg_3,  num_expected_bkg_4,  num_expected_bkg_5, num_expected_bkg_6 );

	}
	void CountingModel::AddChannel(double num_expected_signal, vector<double> num_expected_bkgs){
		AddChannel("",  num_expected_signal,  num_expected_bkgs);
	}
	void CountingModel::AddChannel(vector<double> num_expected_yields, int signal_processes){
		AddChannel("",  num_expected_yields);
	}
	void CountingModel::AddChannel(string channel_name, vector<double> num_expected_yields, int signal_processes){
		if(_debug>=100)cout<<"AddChannel: nsigpro = "<<signal_processes<<endl;
		if(signal_processes<=0)  {cout<<"ERROR: you add a channel with number of signal_processes <=0 "<<endl; exit(0);}
		if(num_expected_yields.size()<=signal_processes){cout<<"ERROR: you add a channel with no background process"<<endl; exit(0);}
		vector<double> num_expected_signals(&num_expected_yields[0], &num_expected_yields[signal_processes]); 
		vector<double> num_expected_bkgs(&num_expected_yields[signal_processes], &num_expected_yields[num_expected_yields.size()]);
		if(_debug>=100)cout<<"channel "<<channel_name<<": tot.size="<<num_expected_yields.size()<<" bkg.size="<<num_expected_bkgs.size()<<endl;
		AddChannel(channel_name, num_expected_signals,  num_expected_bkgs);
	}
	void CountingModel::AddChannel(string channel_name, vector<double> num_expected_signals, vector<double> num_expected_bkgs){
		int signal_processes = num_expected_signals.size();
		int bkg_processes = num_expected_bkgs.size();
		if(signal_processes<=0)  {cout<<"ERROR: you add a channel with number of signal_processes <=0 "<<endl; exit(0);}
		if(bkg_processes<=0)  {cout<<"ERROR: you add a channel with number of bkg_processes <=0 "<<endl; exit(0);}
		v_sigproc.push_back(signal_processes);

		if(channel_name==""){
			char tmp[256];
			sprintf(tmp, "channel_%d", v_channelname.size());
			channel_name==tmp;
			v_channelname.push_back(channel_name);

		}
		else v_channelname.push_back(channel_name);

		double tmp_totbkg = 0;

		vector<double> vsigbkgs; vsigbkgs.clear();
		vector< vector< vector<double> > > vvvuncpar; vvvuncpar.clear();
		vector< vector<int> > vvpdftype; vvpdftype.clear();
		vector< vector<int> > vvidcorrl; vvidcorrl.clear();

		vector< vector<double> > vvunc; vvunc.clear();
		vector<int> vpdftype; vpdftype.clear();
		vector<int> vidcorrl; vidcorrl.clear();

		for(int i=0; i<num_expected_signals.size(); i++){
			vsigbkgs.push_back(num_expected_signals[i]);
			vvvuncpar.push_back(vvunc);
			vvpdftype.push_back(vpdftype);
			vvidcorrl.push_back(vidcorrl);
		}

		for(int i=0; i<num_expected_bkgs.size(); i++){
			vsigbkgs.push_back(num_expected_bkgs[i]);
			vvvuncpar.push_back(vvunc);
			vvpdftype.push_back(vpdftype);
			vvidcorrl.push_back(vidcorrl);
			tmp_totbkg+=num_expected_bkgs[i];
		}

		vv_exp_sigbkgs.push_back(vsigbkgs);
		vv_exp_sigbkgs_scaled.push_back(vsigbkgs);
		v_data.push_back(tmp_totbkg);	
		vvvv_uncpar.push_back(vvvuncpar);
		vvv_pdftype.push_back(vvpdftype);
		vvv_idcorrl.push_back(vvidcorrl);
		if(_debug>=100)cout<<"channel "<<channel_name<<": tot.size="<<signal_processes+bkg_processes<<" bkg.size="<<bkg_processes<<endl;
		vector<string> vproc; vproc.clear();
	       	for(int i=0; i<num_expected_signals.size(); i++) {
			char tmp[256];
			sprintf(tmp, "%d", i-num_expected_signals.size()+1);
			vproc.push_back(tmp);
		}
	       	for(int i=0; i<num_expected_bkgs.size(); i++) {
			char tmp[256];
			sprintf(tmp, "%d", i+1);
			vproc.push_back(tmp);
		}
		vv_procname.push_back(vproc);
	}
	void CountingModel::AddChannel(string channel_name, double num_expected_signal, vector<double> num_expected_bkgs){
		vector<double> num_expected_signals; num_expected_signals.clear();
		num_expected_signals.push_back(num_expected_signal);
		AddChannel(channel_name, num_expected_signals, num_expected_bkgs);
	}	
	void CountingModel::AddUncertainty(int index_channel, int index_sample, double uncertainty_in_relative_fraction, int pdf_type, string uncname ){
		AddUncertainty(index_channel, index_sample, uncertainty_in_relative_fraction, uncertainty_in_relative_fraction, pdf_type, uncname);
	}
	void CountingModel::AddUncertainty(string c, int index_sample, double uncertainty_in_relative_fraction, int pdf_type, string uncname ){
		int index_channel = -1;
		for(int i=0; i<v_channelname.size(); i++){
			if(v_channelname[i]==c) index_channel=i;
		}
		AddUncertainty(index_channel, index_sample, uncertainty_in_relative_fraction, uncertainty_in_relative_fraction, pdf_type, uncname);
	}
	void CountingModel::AddUncertainty(int index_channel, int index_sample, double uncertainty_in_relative_fraction_down, double uncertainty_in_relative_fraction_up, int pdf_type, string uncname ){
		int index_correlation = -1; // numeration starts from 1
		for(int i=0; i<v_uncname.size(); i++){
			if(v_uncname[i]==uncname){
				index_correlation = i+1;	
				break;
			}
		}
		if(index_correlation<0)  {
			index_correlation = v_uncname.size()+1;
			v_uncname.push_back(uncname);
		}
		AddUncertainty(index_channel, index_sample, uncertainty_in_relative_fraction_down, uncertainty_in_relative_fraction_up, pdf_type, index_correlation );

	}
	void CountingModel::AddUncertainty(string c, int index_sample, double uncertainty_in_relative_fraction_down, double uncertainty_in_relative_fraction_up, int pdf_type, string uncname ){
		int index_channel = -1;
		for(int i=0; i<v_channelname.size(); i++){
			if(v_channelname[i]==c) index_channel=i;
		}
		AddUncertainty(index_channel, index_sample, uncertainty_in_relative_fraction_down, uncertainty_in_relative_fraction_up, pdf_type, uncname);
	}
	void CountingModel::AddUncertainty(int index_channel, int index_sample, int npar, double *par, int pdf_type, string uncname ){
		int index_correlation = -1; // numeration starts from 1
		for(int i=0; i<v_uncname.size(); i++){
			if(v_uncname[i]==uncname){
				index_correlation = i+1;	
				break;
			}
		}
		if(index_correlation<0)  {
			index_correlation = v_uncname.size()+1;
			v_uncname.push_back(uncname);
		}
		AddUncertainty(index_channel, index_sample, npar, par, pdf_type, index_correlation );

	}
	void CountingModel::AddUncertainty(string c, int index_sample, int npar, double *par, int pdf_type, string uncname ){
		int index_channel = -1;
		for(int i=0; i<v_channelname.size(); i++){
			if(v_channelname[i]==c) index_channel=i;
		}
		AddUncertainty(index_channel, index_sample, npar, par, pdf_type, uncname);
	}
	void CountingModel::AddUncertainty(int index_channel, int index_sample, double rho, double rho_err, double B, int pdf_type, string uncname ){
		int index_correlation = -1; // numeration starts from 1
		for(int i=0; i<v_uncname.size(); i++){
			if(v_uncname[i]==uncname){
				index_correlation = i+1;	
				break;
			}
		}
		if(index_correlation<0)  {
			index_correlation = v_uncname.size()+1;
			v_uncname.push_back(uncname);
		}
		AddUncertainty(index_channel, index_sample, rho, rho_err, B, pdf_type, index_correlation );
	}
	void CountingModel::AddUncertainty(string c, int index_sample, double rho, double rho_err, double B, int pdf_type, string uncname ){
		int index_channel = -1;
		for(int i=0; i<v_channelname.size(); i++){
			if(v_channelname[i]==c) index_channel=i;
		}
		AddUncertainty(index_channel, index_sample, rho, rho_err, B, pdf_type, uncname);
	}

	// LogNormal and TruncatedGaussian 
	void CountingModel::AddUncertainty(int index_channel, int index_sample, double uncertainty_in_relative_fraction,  int pdf_type, int index_correlation ){
		// sysmetric uncertainties
		AddUncertainty(index_channel, index_sample, uncertainty_in_relative_fraction, uncertainty_in_relative_fraction, pdf_type, index_correlation );
	}
	void CountingModel::AddUncertainty(int index_channel, int index_sample, double uncertainty_in_relative_fraction_down, double uncertainty_in_relative_fraction_up, int pdf_type, int index_correlation ){
		// to deal with asymetric uncertainties
		if( uncertainty_in_relative_fraction_down < 0 or uncertainty_in_relative_fraction_up < 0 ) {
			if(pdf_type==typeTruncatedGaussian) {}; //fine
			if( (uncertainty_in_relative_fraction_down <-1 or uncertainty_in_relative_fraction_up <-1) && pdf_type==typeLogNormal) { cout<<"logNormal type uncertainties can't have kappa < 0, exit"<<endl; exit(0);}; //fine
		} 
		if(pdf_type!=typeLogNormal && pdf_type!= typeTruncatedGaussian && pdf_type!=typeGamma ) {
			cout<<"Error: Currently only implemented LogNormal, Gamma and TruncatedGaussian. Your input "<<pdf_type<<" haven't been implemented, exit"<<endl;
			exit(0);
		}
		if(index_correlation <= 0) { 
			cout<<"Error: index_correlation < 0 "<<endl;
			exit(0);
		}
		vector<double> vunc; vunc.clear(); 
		vunc.push_back(uncertainty_in_relative_fraction_down);
		vunc.push_back(uncertainty_in_relative_fraction_up);
		vvvv_uncpar.at(index_channel).at(index_sample).push_back(vunc);
		vvv_pdftype.at(index_channel).at(index_sample).push_back(pdf_type);
		vvv_idcorrl.at(index_channel).at(index_sample).push_back(index_correlation);
		//ConfigUncertaintyPdfs();
	}
	void CountingModel::AddUncertainty(int index_channel, int index_sample, int npar, double *par, int pdf_type, int index_correlation ){
		if(pdf_type!=typeShapeGaussianQuadraticMorph  && pdf_type!=typeShapeGaussianLinearMorph ) {
			cout<<"Error: AddUncertainty for typeShapeGaussianLinearMorph and typeShapeGaussianQuadraticMorph only. exit"<<endl;
			exit(0);
		}
		if(index_correlation <= 0) { 
			cout<<"Error: index_correlation < 0 "<<endl;
			exit(0);
		}
		vector<double> vunc; vunc.clear(); 
		for(int i=0; i<npar; i++) vunc.push_back(par[i]);
		vvvv_uncpar.at(index_channel).at(index_sample).push_back(vunc);
		vvv_pdftype.at(index_channel).at(index_sample).push_back(pdf_type);
		vvv_idcorrl.at(index_channel).at(index_sample).push_back(index_correlation);
		//ConfigUncertaintyPdfs();
	}

	// From SideBand
	// when B is large, it's close to Gaussian. then use Gaussian, don't use Gamma function
	void CountingModel::AddUncertainty(int index_channel, int index_sample, double rho, double rho_err, double B, int pdf_type, int index_correlation ){
		if(B<0){
			cout<<"Gamma PDF, B can't be less 0"<<endl;
			exit(0);
			//return;
		}
		if(index_correlation <= 0) {
			cout<<"Error: Adding index_correlation <=0, exit "<<endl;
			exit(0);
		}
		if(pdf_type!=typeControlSampleInferredLogNormal && pdf_type!=typeGamma) {
			cout<<"Error: your pdf_type is wrong "<<endl;
			exit(0);
		}
		vector<double> vunc; vunc.clear();
		vunc.push_back(rho);  // if rho<0, it means this gamma term is for multiplicative gamma function...
		vunc.push_back(rho_err);
		vunc.push_back(B);
		vvvv_uncpar.at(index_channel).at(index_sample).push_back(vunc);
		vvv_pdftype.at(index_channel).at(index_sample).push_back(pdf_type);
		vvv_idcorrl.at(index_channel).at(index_sample).push_back(index_correlation);
		//ConfigUncertaintyPdfs();
	}
	void CountingModel::AddObservedData(int index_channel, double num_data){
		v_data.at(index_channel)=num_data;
	}
	void CountingModel::AddObservedData(string c, double num_data){
		int index_channel = -1;
		for(int i=0; i<v_channelname.size(); i++){
			if(v_channelname[i]==c) index_channel=i;
		}
		v_data.at(index_channel)=num_data;
	}

	void CountingModel::SetProcessNames(int index_channel, vector<string> vproc){
		vv_procname.at(index_channel)=vproc;
	}
	void CountingModel::SetProcessNames(string c, vector<string> vproc){
		int index_channel = -1;
		for(int i=0; i<v_channelname.size(); i++){
			if(v_channelname[i]==c) index_channel=i;
		}
		vv_procname.at(index_channel)=vproc;
	}
	void CountingModel::ConfigUncertaintyPdfs(){
		v_TruncatedGaussian_maxUnc.clear();
		v_pdftype.clear();
		max_uncorrelation = 0;
		for(int ch=0; ch<vvv_idcorrl.size(); ch++){
			for(int isam=0; isam<vvv_idcorrl.at(ch).size(); isam++){
				for(int iunc=0; iunc<vvv_idcorrl.at(ch).at(isam).size(); iunc++){
					int indexcorrl = vvv_idcorrl.at(ch).at(isam).at(iunc);
					if(max_uncorrelation<indexcorrl) max_uncorrelation=indexcorrl;
				}
			}
		}
		for(int ch=0; ch<vvv_pdfs_idcorrl.size(); ch++){
			for(int isam=0; isam<vvv_pdfs_idcorrl.at(ch).size(); isam++){
				for(int iunc=0; iunc<vvv_pdfs_idcorrl.at(ch).at(isam).size(); iunc++){
					int indexcorrl = vvv_pdfs_idcorrl.at(ch).at(isam).at(iunc);
					if(max_uncorrelation<indexcorrl) max_uncorrelation=indexcorrl;
				}
			}
		}
		for(int i=0; i<v_pdfs_floatParamsIndcorr.size(); i++){
			if(max_uncorrelation < v_pdfs_floatParamsIndcorr[i]) max_uncorrelation=v_pdfs_floatParamsIndcorr[i];
		}

		for(int i=0; i<=max_uncorrelation; i++){
			v_TruncatedGaussian_maxUnc.push_back(-1);
			v_pdftype.push_back(-1);	
			v_GammaN.push_back(-1);
			double tmpmax=-1;
			for(int ch=0; ch<vvv_idcorrl.size(); ch++){
				for(int isam=0; isam<vvv_idcorrl.at(ch).size(); isam++){
					for(int iunc=0; iunc<vvv_idcorrl.at(ch).at(isam).size(); iunc++){
						int indexcorrl = vvv_idcorrl.at(ch).at(isam).at(iunc);
						if(indexcorrl==i && v_pdftype.back()<0 ){
							v_pdftype.back()=vvv_pdftype.at(ch).at(isam).at(iunc);
						}
						if(indexcorrl==i && vvv_pdftype.at(ch).at(isam).at(iunc)== typeTruncatedGaussian ){
							if(tmpmax< fabs(vvvv_uncpar.at(ch).at(isam).at(iunc).at(0)) ) tmpmax=fabs(vvvv_uncpar.at(ch).at(isam).at(iunc).at(0));	
							if(tmpmax< fabs(vvvv_uncpar.at(ch).at(isam).at(iunc).at(1)) ) tmpmax=fabs(vvvv_uncpar.at(ch).at(isam).at(iunc).at(1));	
						} 
						if(indexcorrl==i && vvv_pdftype.at(ch).at(isam).at(iunc)== typeGamma){
							if(vvvv_uncpar.at(ch).at(isam).at(iunc).at(0)>0 && 
									fabs(vvvv_uncpar.at(ch).at(isam).at(iunc).at(0)*vvvv_uncpar.at(ch).at(isam).at(iunc).at(2) - vv_exp_sigbkgs.at(ch).at(isam)) / vvvv_uncpar.at(ch).at(isam).at(iunc).at(2)/vvvv_uncpar.at(ch).at(isam).at(iunc).at(0)>0.2
							  ) {
								cout<<"channel "<<ch<<"th, "<<v_channelname[ch]<<": process "<<isam<<" using gamma pdf, but rho*B!=b, please check"<<endl; 
								cout<< "rho="<<vvvv_uncpar[ch][isam][iunc][0]<<"  B="<<vvvv_uncpar[ch][isam][iunc][2]<<" b="<<vv_exp_sigbkgs[ch][isam]<<endl;
								exit(0);
							}	
							v_GammaN.back()=vvvv_uncpar.at(ch).at(isam).at(iunc).at(2)+1;	

							/*  from Andrey:
							    The typical convention between the number of observed events in the
							    control region B and pdf(b) is what we wrote in the note. E.g. see

							    -) Bob Cosins's note http://arxiv.org/abs/physics/0702156v4

							    -) Stat commitee recomendation at
http://www.physics.ucla.edu/~cousins/stats/cousins_lognormal_prior.pdf

On the wikipidia page, they derive the pdf starting from a scale-invariant
prior, which I think is a uniform prior for log(B), which is 1/B for B.

If one starts from the uniform prior for B, the answer would be what we
use in the note and in the references above. Such convention makes much
more sense. E.g., the most probable value for b is B*r.
Also, Z_{\Gamma} with the gamma distribution from our note is the same
as the pure frequentist construct Z_{Bi}.

So the bottom line is, let's stick to what we wrote in the note.
If we need to change it later, it will be easy to do.
*/
							//if(v_pdftype.back()==typeGamma)v_GammaN.back()=vvvv_uncpar.at(ch).at(isam).at(iunc).at(2);	
							//						if(v_pdftype.back()==typeGamma)v_GammaN.back()=vvvv_uncpar.at(ch).at(isam).at(iunc).at(2)+1;	
						}
						if(indexcorrl==i && v_pdftype.back()>0 ){
							if( v_pdftype.back()!=vvv_pdftype.at(ch).at(isam).at(iunc) ){
								cout<<" Error:  two uncertainties with 100% correlation must be with same pdftype, exit "<<endl;
								cout<<" Independant unc "<<indexcorrl<<"th, name "<<v_uncname[indexcorrl];
								cout<<" should be "<<v_pdftype.back()<<", but "<<vvv_pdftype.at(ch).at(isam).at(iunc)<<" for ch"<<ch<<"("<<v_channelname[ch]<<")"
									<<" isam"<<isam<<"("<<vv_procname[ch][isam]<<")"
									<<" iunc"<<iunc<<"("<<v_uncname[indexcorrl]<<")"<<endl;
								cout<<" The conflict unc at above process is for "<<vvv_idcorrl[ch][isam][iunc]
									<<"th source, name "<<v_uncname[vvv_idcorrl[ch][isam][iunc]]<<endl;
								exit(0);
							}
						}
					}
				}
			}
			for(int ch=0; ch<vvv_pdfs_idcorrl.size(); ch++){
				for(int isam=0; isam<vvv_pdfs_idcorrl.at(ch).size(); isam++){
					for(int iunc=0; iunc<vvv_pdfs_idcorrl.at(ch).at(isam).size(); iunc++){
						int indexcorrl = vvv_pdfs_idcorrl.at(ch).at(isam).at(iunc);
						if(indexcorrl==i && v_pdftype.back()<0 ){
							v_pdftype.back()=vvv_pdfs_pdftype.at(ch).at(isam).at(iunc);
						}
						if(indexcorrl==i && vvv_pdfs_pdftype.at(ch).at(isam).at(iunc)== typeTruncatedGaussian ){
							if(tmpmax< fabs(vvv_pdfs_normvariation.at(ch).at(isam).at(iunc).at(0)) ) tmpmax=fabs(vvv_pdfs_normvariation.at(ch).at(isam).at(iunc).at(0));	
							if(tmpmax< fabs(vvv_pdfs_normvariation.at(ch).at(isam).at(iunc).at(1)) ) tmpmax=fabs(vvv_pdfs_normvariation.at(ch).at(isam).at(iunc).at(1));	
						} 
						if(indexcorrl==i && vvv_pdfs_pdftype.at(ch).at(isam).at(iunc)== typeGamma){
							if(vvv_pdfs_normvariation.at(ch).at(isam).at(iunc).at(0)>0 && 
									fabs(vvv_pdfs_normvariation.at(ch).at(isam).at(iunc).at(0)*vvv_pdfs_normvariation.at(ch).at(isam).at(iunc).at(2) - vv_pdfs_norm.at(ch).at(isam)) / vvv_pdfs_normvariation.at(ch).at(isam).at(iunc).at(2)/vvv_pdfs_normvariation.at(ch).at(isam).at(iunc).at(0)>0.2
							  ) {
								cout<<"Shape channel "<<ch<<"th, "<<v_pdfs_channelname[ch]<<": process "<<isam<<" using gamma pdf, but rho*B!=b, please check"<<endl; 
								cout<< "rho="<<vvv_pdfs_normvariation[ch][isam][iunc][0]<<"  B="<<vvv_pdfs_normvariation[ch][isam][iunc][2]<<" b="<<vv_exp_sigbkgs[ch][isam]<<endl;
								exit(0);
							}	
							v_GammaN.back()=vvv_pdfs_normvariation.at(ch).at(isam).at(iunc).at(2)+1;	

						}
						if(indexcorrl==i && v_pdftype.back()>0 ){
							if( v_pdftype.back()!=vvv_pdfs_pdftype.at(ch).at(isam).at(iunc) ){
								cout<<" Error:  two uncertainties with 100% correlation must be with same pdftype, exit "<<endl;
								cout<<" Independant unc "<<indexcorrl<<"th, name "<<v_uncname[indexcorrl];
								cout<<" should be "<<v_pdftype.back()<<", but "<<vvv_pdfs_pdftype.at(ch).at(isam).at(iunc)<<" for shape ch"<<ch<<"("<<v_pdfs_channelname[ch]<<")"
									<<" isam"<<isam<<"("<<vv_pdfs[ch][isam]<<")"
									<<" iunc"<<iunc<<"("<<v_uncname[indexcorrl]<<")"<<endl;
								cout<<" The conflict unc at above process is for "<<vvv_pdfs_idcorrl[ch][isam][iunc]
									<<"th source, name "<<v_uncname[vvv_idcorrl[ch][isam][iunc]]<<endl;
								exit(0);
							}
						}
					}
				}
			}
			if(tmpmax>0){
				v_TruncatedGaussian_maxUnc.back()=tmpmax;
			} 
		}
		for(int i=0; i<v_pdfs_floatParamsName.size(); i++){
			v_pdftype[v_pdfs_floatParamsIndcorr[i]] = typeBifurcatedGaussian;
		}
		MakeListOfShapeUncertainties();
	}	
	void CountingModel::MakeListOfShapeUncertainties(){
		vvv_shapeuncindex.clear();
		for(int ch=0; ch<vvv_idcorrl.size(); ch++){
			vector< vector<int> > vvp; vvp.clear();
			for(int isam=0; isam<vvv_idcorrl.at(ch).size(); isam++){
				vector<int> vshape; vshape.clear();
				for(int iunc=0; iunc<vvv_idcorrl.at(ch).at(isam).size(); iunc++){
					int indexcorrl = vvv_idcorrl.at(ch).at(isam).at(iunc);
					if(v_pdftype[indexcorrl]==typeShapeGaussianQuadraticMorph or v_pdftype[indexcorrl]==typeShapeGaussianLinearMorph)
						vshape.push_back(iunc);
				}
				vvp.push_back(vshape);
			}
			vvv_shapeuncindex.push_back(vvp);
		}

	}
	VChannelVSample CountingModel::FluctuatedNumbers(double *par){
		if(_rdm==NULL) {cout<<"Model random gen engine not set yet, exit "<<endl; exit(0);}
		if(!b_systematics) {
			if(bHasParametricShape){
				vv_pdfs_params_varied = vv_pdfs_params;
				vv_pdfs_norm_varied = vv_pdfs_norm_scaled;
				for(int ch=0; ch<vv_pdfs.size(); ch++){
					int nsigproc = v_pdfs_sigproc[ch];
					for(int isam=0; isam<vv_pdfs[ch].size(); isam++){
						_w_varied->var(vv_pdfs_normNAME[ch][isam])->setVal(vv_pdfs_norm_varied[ch][isam]);
					}
				}
			}
			return vv_exp_sigbkgs_scaled;
		}

		double tmp ; 
		vector<double> vrdm; vrdm.clear();
		v_pdfs_floatParamsVaried.clear();
		if(par==0){
			for(int i=0; i<v_pdftype.size(); i++){
				if(_debug>=100)cout<<" vpdftype: "<<i<<"th --> "<<v_pdftype[i]<<endl;
				vrdm.push_back(-999);
				switch (v_pdftype[i]){
					case typeLogNormal:
						vrdm.back()=_rdm->Gaus();
						break;
					case typeShapeGaussianLinearMorph:
					case typeShapeGaussianQuadraticMorph:
						/*	
							tmp = -5;
							while(fabs(tmp)>4){
							tmp=_rdm->Gaus();
							}
							vrdm.back()=tmp;
							break;
							*/	
						vrdm.back()=_rdm->Gaus();
						break;

					case typeTruncatedGaussian:
						//      another way is to throw normal gaus random number and regenerate if x<-1, it's more transparent
						tmp = -2;
						while(tmp<-1){
							tmp=_rdm->Gaus(0, v_TruncatedGaussian_maxUnc[i]);
						}
						vrdm.back()=tmp;
						break;

					case typeGamma:
						//if(_debug)cout<<" i = "<<i<<"  v_GammaN[i]="<<v_GammaN[i]<<endl;
						vrdm.back()=_rdm->Gamma(v_GammaN[i]);
						//if(_debug) cout<<"done for random gamma"<<endl;
						break;
					case typeBifurcatedGaussian:
						{
							if(_debug>=100) cout<<" generating new value for parameter "<<v_uncname[i-1]<<endl;
							RooDataSet *tmpRDS =  _w_varied->pdf(TString::Format("%s_bfg",v_uncname[i-1].c_str()))
								->generate(RooArgSet(*_w_varied->var(TString::Format("%s_x", v_uncname[i-1].c_str()))), 1);
							vrdm.back()= dynamic_cast<RooRealVar*> ( tmpRDS->get(0)->first() )->getVal();
							delete tmpRDS;
							if(_debug>=100) cout<<"  "<<vrdm.back()<<endl;
							v_pdfs_floatParamsVaried.push_back(vrdm.back());
							break;
						}
					case typeControlSampleInferredLogNormal:
						//dummy
						cout<<"Error: We haven't implemented the pdf of typeControlSampleInferredLogNormal"<<endl;
						exit(0);
						break;
					default:
						break;
						//dummy
						//cout<<"Error: Unknown pdf_type "<<v_pdftype[i]<<endl;
						//exit(0);
				}	
				//if(_debug) cout<<"done for random gen "<<i<<endl;
			}
		}else{
			vrdm.push_back(0);
			for(int i=1; i<v_pdftype.size(); i++){
				vrdm.push_back(par[i]);
				if(_debug>=100)cout<<" index "<<i<<": "<<par[i]<<endl;
			}
		}		
		//if(_debug) cout<<"done for random gen"<<endl;

		VChannelVSample vv = vv_exp_sigbkgs_scaled;
		int indexcorrl, pdftype, isam, iunc;
		vector<int> shapeuncs;
		double ran = 0, h;
		//if(_debug) cout<<"vvv_idcorrl.size="<<vvv_idcorrl.size()<<endl;
		int nsigproc = 1;
		double * uncpars;
		double tmprand;
		bool added = false;
		double norminal = 0;
		double normalization = 0;
		for(int ch=0; ch<vvv_idcorrl.size(); ch++){
			nsigproc = v_sigproc[ch];
			for(isam=0; isam<vvv_idcorrl[ch].size(); isam++){
				if(bMoveUpShapeUncertainties){
					shapeuncs = vvv_shapeuncindex[ch][isam];
					h=0;
					added = false;
					for(int i = 0; i<shapeuncs.size(); i++){
						indexcorrl = vvv_idcorrl[ch][isam][shapeuncs[i]];
						pdftype = vvv_pdftype[ch][isam][shapeuncs[i]];
						ran = vrdm[indexcorrl];
						uncpars  = &(vvvv_uncpar[ch][isam][shapeuncs[i]][0]);
						switch (pdftype){
							case typeShapeGaussianLinearMorph:
								if(*(uncpars+7) == 1.){
									tmprand = ran; 
									ran*= (*(uncpars+6));
									if(!added) {h+=*(uncpars+2); added=true; norminal = h; normalization = *(uncpars+3); }
									//h += max(-ran, 0.) * (*uncpars) + max(ran, 0.) * (*(uncpars+1)) - fabs(ran)*(*(uncpars+2)) ; // uncer params:  down, up, norminal, normlization_of_main_histogram,  uncertainty_down_onNorm, uncertainty_up_onNorm, scale_down_of_gaussian,  siglebin_or_binned
									h += ( ran * (ran>0? *(uncpars+1)-*(uncpars+2): *(uncpars+2)-*uncpars) );
									ran = tmprand;
								}

								break;
							case typeShapeGaussianQuadraticMorph:
								if(*(uncpars+7) == 1.){
									tmprand = ran; 
									ran*= (*(uncpars+6));
									if(!added) {h+=*(uncpars+2); added=true; norminal = h;  normalization = *(uncpars+3); }
									if(fabs(ran)<1){
										h += ran * (ran-1)/2. * (*uncpars) + ran * (ran+1)/2. * (*(uncpars+1)) - ran*ran*(*(uncpars+2)) ; 
									} else { 
										//h += max(-ran, 0.) * (*uncpars) + max(ran, 0.) * (*(uncpars+1)) - fabs(ran)*(*(uncpars+2)) ; 
										h += ( ran * (ran>0? *(uncpars+1)-*(uncpars+2): *(uncpars+2)-*uncpars) );
									}
									ran = tmprand;
								}
								break;
							default:
								break;
						}
					}

					if(added){
						if(h<=0) h=10e-9;
						/*
						   if( norminal !=0 && vv[ch][isam]!=0) {
						   vv[ch][isam]*=h/norminal;	
						   }else if(vv[ch][isam]==0) vv[ch][isam] = h*normalization;
						   else { ;}
						   */

						if(_debug>=100)	cout<<"c="<<ch<<" s="<<isam<<" normalization = "<<normalization<<" ran="<<ran
							<<" h="<<h<<" h*normalization="<<h*normalization
								<<" bs="<<vv[ch][isam]<<" norminal="<<norminal<<" h/norminal="<<h/norminal
								<<" bs*=h/norminal = "<<(vv[ch][isam]*(h/norminal)) <<endl;	
						if( norminal !=0) 
							vv[ch][isam]*=h/norminal;	
						else
							vv[ch][isam]=h*normalization;
						//here no need to take again the r into account....   if norminal = 0
					}
				}

				for(iunc=0; iunc<vvv_idcorrl[ch][isam].size(); iunc++){
					//if(_debug) cout<<ch<<" "<<isam<<" "<<iunc<<" "<<endl;
					indexcorrl = vvv_idcorrl[ch][isam][iunc];
					pdftype = vvv_pdftype[ch][isam][iunc];
					ran = vrdm[indexcorrl];
					uncpars  = &(vvvv_uncpar[ch][isam][iunc][0]);
					switch (pdftype){
						case typeLogNormal : 
							//vv[ch][isam]*=pow( (1+ vvvv_uncpar[ch][isam][iunc][ (ran>0?1:0) ]), ran ); // down/up
							vv[ch][isam]*=pow( (1+ (ran>0? *(uncpars+1):*uncpars) ) , ran );
							break;

						case typeTruncatedGaussian :
							vv[ch][isam]*=( 1+vvvv_uncpar[ch][isam][iunc][(ran>0?1:0)] / v_TruncatedGaussian_maxUnc[indexcorrl] * ran );		
							break;

						case typeGamma :
							if(vvvv_uncpar[ch][isam][iunc][0]>0){
								tmp = vv_exp_sigbkgs_scaled[ch][isam];
								if(isam<nsigproc){
									if(tmp!=0) {vv[ch][isam] /=tmp; vv[ch][isam]*=(ran * vvvv_uncpar[ch][isam][iunc][0] * _common_signal_strength );}
									else vv[ch][isam] = ran * vvvv_uncpar[ch][isam][iunc][0] * _common_signal_strength ; // Gamma
								}else{
									if(tmp!=0) {vv[ch][isam] /=tmp; vv[ch][isam]*=(ran * vvvv_uncpar[ch][isam][iunc][0] );}
									else vv[ch][isam] = ran * vvvv_uncpar[ch][isam][iunc][0]; // Gamma
								}
							}else{ // if rho<0,   then this is multiplicative gamma function ....
								vv[ch][isam] *= (ran/v_GammaN[indexcorrl]);
							}
							break;

						case typeControlSampleInferredLogNormal :
							// FIXME
							break;
						case typeShapeGaussianLinearMorph:
							if(!bMoveUpShapeUncertainties){
								if(*(uncpars+7) == 1.){
									tmprand = ran; 
									ran*= (*(uncpars+6));
									h = *(uncpars+2) + max(-ran, 0.) * (*uncpars) + max(ran, 0.) * (*(uncpars+1)) - fabs(ran)*(*(uncpars+2)) ; // uncer params:  down, up, norminal, normlization_of_main_histogram,  uncertainty_down_onNorm, uncertainty_up_onNorm, scale_down_of_gaussian,  siglebin_or_binned
									//if(h<0) h=0;
									if(h<=0) h=10e-9;
									//if(vv_exp_sigbkgs_scaled[ch][isam]!=0 && vv[ch][isam]!=0) 
									if( *(uncpars+2)!=0 && vv[ch][isam]!=0) {
										vv[ch][isam]*=h/(*(uncpars+2));	
									}else if(vv[ch][isam]==0) vv[ch][isam] = (*(uncpars+3))*h;
									else { ;}
									ran = tmprand;
								}
							}
							if(_debug>=100)cout<< "main =" << *(uncpars+2) << " up="<< *(uncpars+1) << " down="<< *uncpars << " ran="<<ran<< " --> h="<<h<<endl;
							//vv[ch][isam]*=pow( (ran>0? *(uncpars+5):*(uncpars+4)) , ran*(*(uncpars+6)) );
							vv[ch][isam]*=pow( (ran>0? *(uncpars+5):*(uncpars+4)) , ran>0?ran*(*(uncpars+6)): -ran*(*(uncpars+6)));
							if(_debug>=100)cout<<ch<<" "<<isam<<" "<<iunc<<" : "<<vv[ch][isam]<<endl;

							break;
						case typeShapeGaussianQuadraticMorph:
							//vv[ch][isam]*=pow( (ran>0? *(uncpars+5):*(uncpars+4) ) , ran*(*(uncpars+6)) );
							vv[ch][isam]*=pow( (ran>0? *(uncpars+5):*(uncpars+4)) , ran>0?ran*(*(uncpars+6)): -ran*(*(uncpars+6)));
							if(!bMoveUpShapeUncertainties){
								if(*(uncpars+7) == 1.){
									tmprand = ran; 
									ran*= (*(uncpars+6));
									if(fabs(ran)<1)
										h = *(uncpars+2) + ran * (ran-1)/2. * (*uncpars) + ran * (ran+1)/2. * (*(uncpars+1)) - ran*ran*(*(uncpars+2)) ; // uncer params:  down, up, norminal, normlization_of_main_histogram,  uncertainty_down_onNorm, uncertainty_up_onNorm, scale_down_of_gaussian, siglebin_or_binned
									else 
										h = *(uncpars+2) + max(-ran, 0.) * (*uncpars) + max(ran, 0.) * (*(uncpars+1)) - fabs(ran)*(*(uncpars+2)) ; // uncer params:  down, up, norminal, normlization_of_main_histogram,  uncertainty_down_onNorm, uncertainty_up_onNorm,  scale_down_of_gaussian, siglebin_or_binned
									//if(h<0) h=0;
									if(h<=0) h=10e-9;
									if( *(uncpars+2)!=0 && vv[ch][isam]!=0) {
										vv[ch][isam]*=h/(*(uncpars+2));	
									}else if(vv[ch][isam]==0) vv[ch][isam] = (*(uncpars+3))*h;
									else { ;}
									ran = tmprand;
								}
							}
							if(_debug>=100)cout<< "main =" << *(uncpars+2) << " up="<< *(uncpars+1) << " down="<< *uncpars << " ran="<<ran<< " --> h="<<h<<endl;


							break;
						default:
							break;
					}

				}
			}
		}

		if(bHasParametricShape){
			vv_pdfs_params_varied = vv_pdfs_params;
			vv_pdfs_norm_varied = vv_pdfs_norm_scaled;
			for(int ch=0; ch<vv_pdfs.size(); ch++){
				nsigproc = v_pdfs_sigproc[ch];
				for(isam=0; isam<vv_pdfs[ch].size(); isam++){
					for(iunc=0; iunc<vvv_pdfs_idcorrl[ch][isam].size(); iunc++){
						indexcorrl = vvv_pdfs_idcorrl[ch][isam][iunc];
						pdftype = vvv_pdfs_pdftype[ch][isam][iunc];
						ran = vrdm[indexcorrl];
						switch (pdftype){
							default:
								if(_debug>=100) cout<<" in FluctuatedNumbers, bHasParametricShape: ch["<<ch<<"] isam["<<isam<<"] iunc["<<iunc<<"]"<<" pdftype="<<pdftype<<endl;
								/*
								   for(int i=0; i<vv_pdfs_npar[ch][isam]; i++){
								   vv_pdfs_params_varied[ch][isam][i] += 
								   ( ran * (ran>0? vvv_pdfs_params_up[ch][isam][iunc][i] -vv_pdfs_params[ch][isam][i] 
								   : vv_pdfs_params[ch][isam][i] -vvv_pdfs_params_down[ch][isam][iunc][i] ) );
								   }
								   */		
								vv_pdfs_norm_varied[ch][isam] *= pow( (1+ (ran>0? vvv_pdfs_normvariation[ch][isam][iunc][1]: vvv_pdfs_normvariation[ch][isam][iunc][0])) , ran );
								if(_debug>=100) cout<<" norm varied after this unc = "<<vv_pdfs_norm_varied[ch][isam]<<endl;
								break;

						}
					}
					if(_debug>=100) cout<<" norm varied after all unc = "<<vv_pdfs_norm_varied[ch][isam]<<endl;
					_w_varied->var(vv_pdfs_normNAME[ch][isam])->setVal(vv_pdfs_norm_varied[ch][isam]);
				}
			}
			for(int i=0; i<v_pdfs_floatParamsName.size(); i++){
				if(_debug>=100) cout<<" setting new value for parameter "<<v_pdfs_floatParamsName[i]<<": "<<vrdm[v_pdfs_floatParamsIndcorr[i]]<<endl;
				_w_varied->var(v_pdfs_floatParamsName[i].c_str())->setVal(vrdm[v_pdfs_floatParamsIndcorr[i]]);
			}
			if(_debug>=100) {
				cout<<"FluctuatedNumbers, varied workspace: "<<endl;
				//_w_varied->Print("V"); // cause crash   ..... at some point
			}
		}

		return vv;
	}	
	VIChannel CountingModel::GetToyData_H0(){
		// background only hypothesis
		VChannelVSample vv = FluctuatedNumbers();

		VIChannel v; v.clear();
		double tmp;
		for(int ch=0; ch<vv.size(); ch++){
			tmp=0;
			for(int isam=v_sigproc[ch]; isam<vv[ch].size(); isam++){
				//start from 1,  don't add signal
				tmp+=vv[ch][isam];	
			}
			v.push_back(_rdm->Poisson(tmp));
		}

		if(bHasParametricShape){
			for(int ch=0; ch<v_pdfs_roodataset_toy.size(); ch++){
				delete v_pdfs_roodataset_toy[ch];
			}
			v_pdfs_roodataset_toy.clear();
			for(int ch=0; ch<vv_pdfs.size(); ch++){
				vector<double> vtmp;vtmp.clear();
				if(_debug>=100)cout<<" in GetToyData_H0, bHasParametricShape: ch = "<<ch<<endl;
				v_pdfs_roodataset_toy.push_back( _w_varied->pdf(v_pdfs_b[ch])->generate(RooArgSet(*_w->var(v_pdfs_observables[ch])), Extended()));

				if(_debug>=1000)v_pdfs_roodataset_toy[ch]->Print("V");

				if(_debug>=1000){
					TString s = "H0_"; s+= (int)(10000000*_rdm->Rndm()); s+=".root";
					TFile f(s, "RECREATE") ;
					TH1F h("H0","H0", 100, _w->var(v_pdfs_observables[ch])->getMin(), _w->var(v_pdfs_observables[ch])->getMax() );
					v_pdfs_roodataset_toy[ch]->fillHistogram(&h, RooArgList(*_w->var(v_pdfs_observables[ch])));
					f.WriteTObject(&h);
					f.Close();
				}

				if(_debug>=10)cout<<"H0, number of events generated for channel "<<ch<<": "<<v_pdfs_roodataset_toy[ch]->sumEntries()<<endl;
				for(int i=0; i<v_pdfs_roodataset_toy[ch]->sumEntries(); i++){
					RooRealVar *r = dynamic_cast<RooRealVar*>(v_pdfs_roodataset_toy[ch]->get(i)->first());
					vtmp.push_back( r->getVal() );
				}
				vv_pdfs_data_toy.push_back(vtmp);
			}
		}	

		return v;
	}
	VIChannel CountingModel::GetToyData_H1(){
		// alternative hypothesis
		VChannelVSample vv = FluctuatedNumbers();

		if(_debug>=100)cout<<" in GetToyData_H1 ......"<<endl;

		VIChannel v; v.clear();
		double tmp;
		for(int ch=0; ch<vv.size(); ch++){
			tmp=0;
			for(int isam=0; isam<vv[ch].size(); isam++){
				//start from 0, add signal
				tmp+=vv[ch][isam];	
			}
			v.push_back(_rdm->Poisson(tmp));
		}

		if(_debug>=100)cout<<" in GetToyData_H1 done for counting channels"<<endl;
		if(bHasParametricShape){
			if(_debug>=100)cout<<" in GetToyData_H1 start to destroy old toy"<<endl;
			for(int ch=0; ch<v_pdfs_roodataset_toy.size(); ch++){
				if(v_pdfs_roodataset_toy[ch])	delete v_pdfs_roodataset_toy[ch];
			}
			v_pdfs_roodataset_toy.clear();
			if(_debug>=100)cout<<" in GetToyData_H1 old toy destroyed"<<endl;
			for(int ch=0; ch<vv_pdfs.size(); ch++){
				vector<double> vtmp;vtmp.clear();

				//FIXME   we need release memory of used toys,  but delete them cause some segmentation fault 
				//if(!v_pdfs_roodataset_toy[ch]->IsZombie()) delete v_pdfs_roodataset_toy[ch];
				//v_pdfs_roodataset_toy[ch]->Delete();

				if(_debug>=100)cout<<" in GetToyData_H1, bHasParametricShape: ch = "<<ch<<endl;
				v_pdfs_roodataset_toy.push_back(_w_varied->pdf(v_pdfs_sb[ch])->generate(RooArgSet(*_w->var(v_pdfs_observables[ch])), Extended()));

				if(_debug>=100)v_pdfs_roodataset_toy[ch]->Print("V");

				if(_debug>=1000){
					TString s = "H1_"; s+= (int)(10000000*_rdm->Rndm()); s+=".root";
					TFile f(s, "RECREATE") ;
					TH1F h("H1","H1", 100, _w->var(v_pdfs_observables[ch])->getMin(), _w->var(v_pdfs_observables[ch])->getMax() );
					v_pdfs_roodataset_toy[ch]->fillHistogram(&h, RooArgList(*_w->var(v_pdfs_observables[ch])));
					h.Write();
					f.Write();
					f.Close();
				}

				if(_debug>=10) cout<<"H1, number of events generated for channel "<<ch<<": "<<v_pdfs_roodataset_toy[ch]->sumEntries()<<endl;
				for(int i=0; i<v_pdfs_roodataset_toy[ch]->sumEntries(); i++){
					RooRealVar *r = dynamic_cast<RooRealVar*>(v_pdfs_roodataset_toy[ch]->get(i)->first());
					vtmp.push_back( r->getVal() );
				}
				vv_pdfs_data_toy.push_back(vtmp);
			}
		}	

		return v;
	}
	void CountingModel::Print(int printLevel){

		if(printLevel >= 100) {
			// this print out only work for  nsigproc=1

			for(int i=0; i<v_sigproc.size(); i++){
				if(v_sigproc[i]>1) return;
			}

			//	cout<<"I'm here, dummy"<<endl;
			cout<<"\n\n\t ***********Start Model Printout*************"<<endl;		
			cout<<"  -------- print out, messy version  ---------"<<endl;

			if(_common_signal_strength!=1) cout<<"\t signal rates have been scaled by a factor of "<<_common_signal_strength<<endl;
			cout<<"nevents   ";	
			for(int ch=0; ch<vv_exp_sigbkgs_scaled.size(); ch++){
				for(int isam=0; isam<vv_exp_sigbkgs_scaled[ch].size(); isam++){
					printf("%8.3f", vv_exp_sigbkgs_scaled[ch][isam] );
				}
			}
			if(b_systematics) {
				cout<<endl<<endl;
				cout<<"uncertn   ";	
				for(int ch=0; ch<vvvv_uncpar.size(); ch++){
					for(int isam=0; isam<vvvv_uncpar[ch].size(); isam++){
						if(vvvv_uncpar[ch][isam].size()<=0) continue;
						printf("%8.3f", vvvv_uncpar[ch][isam][0][0]);
					}
				}

				cout<<endl<<endl;
				cout<<"correla   ";	
				for(int ch=0; ch<vvv_idcorrl.size(); ch++){
					for(int isam=0; isam<vvv_idcorrl[ch].size(); isam++){
						if(vvv_idcorrl[ch][isam].size()<=0) continue;
						printf("%8d", vvv_idcorrl[ch][isam][0]);
					}
				}
				cout<<endl<<endl;
				cout<<"pdftype   ";	
				for(int ch=0; ch<vvv_pdftype.size(); ch++){
					for(int isam=0; isam<vvv_pdftype[ch].size(); isam++){
						if(vvv_pdftype[ch][isam].size()<=0) continue;
						printf("%8d", vvv_pdftype[ch][isam][0]);
					}
				}
				cout<<endl<<endl;
			}
			//    rotate the table,  simplified version
			cout<<" \n\t -------- rotate print out, simplified version  ---------"<<endl;
			if(_common_signal_strength!=1) cout<<"\t signal rates have been scaled by a factor of "<<_common_signal_strength<<endl;
			if(0 && b_systematics) {
				printf("%8s%8s%8s%8s%8s%8s%8s%8s%8s%8s\n", "channel", "signal", "es", "bkg", "eb", "pdf_es", "pdf_eb", "id_cor", "id_cor", "data");
				for(int ch=0; ch<vv_exp_sigbkgs_scaled.size(); ch++) {
					double tmp_totbkg=0;
					for(int ns=1; ns<vv_exp_sigbkgs_scaled[ch].size(); ns++) {
						tmp_totbkg+=vv_exp_sigbkgs_scaled[ch][ns];	
					}
					printf("%8d%8.3f%8.2f%8.3f%8.2f%8d%8d%8d%8d%8.1f\n",
							ch,
							vv_exp_sigbkgs_scaled[ch][0],
							vvvv_uncpar[ch][0][0][0],
							//vv_exp_sigbkgs_scaled[ch][1],
							tmp_totbkg,
							vvvv_uncpar[ch][1][0][0],
							vvv_pdftype[ch][0][0],
							vvv_pdftype[ch][1][0],
							vvv_idcorrl[ch][0][0],
							vvv_idcorrl[ch][1][0],
							v_data[ch]	
					      );
				}
			}else{
				printf("%12s%12s%12s%12s\n", "channel", "signal", "bkg", "data");
				for(int ch=0; ch<vv_exp_sigbkgs_scaled.size(); ch++) {
					double tmp_totbkg=0;
					for(int ns=1; ns<vv_exp_sigbkgs_scaled[ch].size(); ns++) {
						tmp_totbkg+=vv_exp_sigbkgs_scaled[ch][ns];	
					}
					printf("%12d%12.6f%12.6f%12.6f\n",
							ch,
							vv_exp_sigbkgs_scaled[ch][0],
							//vv_exp_sigbkgs_scaled[ch][1],
							tmp_totbkg,
							v_data[ch]	
					      );
				}
			}
		}
		cout<<" v_uncname.size= "<<v_uncname.size()<<";  v_pdftype.size= "<<v_pdftype.size()<<endl;
		cout<<" list of independant uncertainties: (index  name  type)"<<endl;
		for(int i=0; i<v_uncname.size(); i++){
			cout<<i+1<<" "<<v_uncname[i]<<"  ";
			if(v_pdftype.size()>0) cout	<< v_pdftype[i+1] ;
			cout<<endl;
		}
		cout<<endl;

		//  more detail print out
		cout<<" -------- more detail print out ---------"<<endl;
		int step=0;
		if(printLevel<=1) step= vv_exp_sigbkgs_scaled.size()/10;
		for(int ch=0; ch<vv_exp_sigbkgs_scaled.size(); ch+=(step+1)) {
			cout<<endl;
			cout<<"c"<<ch<<" ("<<v_channelname[ch]<<"):"<<endl;
			for(int ns=0; ns<vv_exp_sigbkgs_scaled[ch].size(); ns++) {
				if(ns<v_sigproc[ch])cout<<"  signal "<<vv_procname[ch][ns]<<": events "<<vv_exp_sigbkgs_scaled[ch][ns]<<endl;
				else cout<<"  bkg "<<vv_procname[ch][ns]<<": events "<<vv_exp_sigbkgs_scaled[ch][ns]<<endl;
				if(b_systematics){	for(int ierr=0; ierr<vvvv_uncpar[ch][ns].size(); ierr++)
					cout<<" \t unc "<<v_uncname[vvv_idcorrl[ch][ns][ierr] - 1]<<" "<< vvv_pdftype[ch][ns][ierr]<< " "  <<" "<<vvvv_uncpar[ch][ns][ierr][0]<<endl;	
				}
			}
			cout<<"\t\t  observed data events = " << v_data[ch]<<endl;
		}


		fflush(stdout);

		cout<<"\t ***********End Model Printout*************\n"<<endl;		
	}
	bool CountingModel::Check(){
		// check any two uncertainties have same index_correlation  are with same pdf_type
		if( v_data.size()!= vv_exp_sigbkgs.size() || v_data.size()!= vvvv_uncpar.size() ){
			cout<<"channels of data != expected signal or backgrounds, or != vvvv_uncpar.size "<<endl;
			cout<<"data.size="<<int(v_data.size())<<endl;
			cout<<"vv_exp_sigbkgs.size="<<(int)vv_exp_sigbkgs.size()<<endl;
			cout<<"vvvv_uncpar.size="<<(int)vvvv_uncpar.size()<<endl;
			exit(0);
		}	
		if(!_rdm) {
			cout<<"random generator not set yet, exit"<<endl;
		}
		for(int ch=0; ch<v_data.size(); ch++){
			if(v_data[ch]<0) {
				cout<<"\t data in channel "<<ch<<" < 0 ("<<v_data[ch]<<"), exit"<<endl;
				exit(0);
			}
		}	

		// need more check
		return true;
	}
	void CountingModel::SetSignalScaleFactor(double r){
		if(_debug>=10) cout<<"\n  *** SetSignalScaleFactor r= "<<r<<endl;
		if(!b_AllowNegativeSignalStrength && r<=0 ){
			cout<<"Error: signal strength r <=0"<<endl;
			cout<<"If you want to allow signal strength to be non-positive, please \n *** model->SetAllowNegativeSignalStrength(true)"<<endl;
			return;
		}
		_common_signal_strength=r;
		vv_exp_sigbkgs_scaled = vv_exp_sigbkgs;
		for(int ch=0; ch<vv_exp_sigbkgs_scaled.size(); ch++){
			for(int isam=0; isam<v_sigproc[ch]; isam++){
				vv_exp_sigbkgs_scaled[ch][isam]*=_common_signal_strength;
			}
		}
		vv_pdfs_norm_scaled = vv_pdfs_norm;
		for(int ch=0; ch<vv_pdfs_norm_scaled.size(); ch++){
			for(int isam=0; isam<v_pdfs_sigproc[ch]; isam++){
				vv_pdfs_norm_scaled[ch][isam]*=_common_signal_strength;
				_w->var(vv_pdfs_normNAME[ch][isam])->setVal(vv_pdfs_norm_scaled[ch][isam]);
				_w_varied->var(vv_pdfs_normNAME[ch][isam])->setVal(vv_pdfs_norm_scaled[ch][isam]);
			}
		}

		if(_debug>=1000 and bHasParametricShape){
			TString s = "pdf_r"; s+=r; s+="_"; s+=".root";
			TFile f(s, "RECREATE") ;
			for(int ch=0; ch<vv_pdfs.size(); ch++){
				double xmin = _w->var(v_pdfs_observables[ch])->getMin();
				double xmax = _w->var(v_pdfs_observables[ch])->getMax();
				TString chnm = v_pdfs_channelname[ch];	
				TH1F hs(chnm+"_s","s", 100, xmin, xmax);
				TH1F hb(chnm+"_b","b", 100,  xmin, xmax);
				TH1F hsb(chnm+"_sb","sb", 100,  xmin, xmax);
				RooArgSet vars(*_w->var(v_pdfs_observables[ch]));
				for(int x = 1; x<=100; x++){
					_w->var(v_pdfs_observables[ch])->setVal((xmax-xmin)/100.*(x-1)+xmin); 
					hs.SetBinContent(x, _w->pdf(v_pdfs_s[ch])->getVal(&vars));
					hb.SetBinContent(x, _w->pdf(v_pdfs_b[ch])->getVal(&vars));
					hsb.SetBinContent(x, _w->pdf(v_pdfs_sb[ch])->getVal(&vars));
				}

				hs.Write();
				hb.Write();
				hsb.Write();
			}
			f.Close();
		}

		//if allow signal strength to be non-positive, then please make sure sig+bkgs >=0 in each channel 
		if(b_AllowNegativeSignalStrength){
			for(int ch=0; ch<vv_exp_sigbkgs_scaled.size(); ch++){
				double tot = 0;
				for(int p=0; p<vv_exp_sigbkgs_scaled[ch].size(); p++)	tot+=vv_exp_sigbkgs_scaled[ch][p];
				/*if(tot<0) {
				  cout<<"Error: negative tot yield in channel "<<ch+1<<endl;
				  cout<<"Please SetAllowNegativeSignalStrength(false)"<<endl;
				  cout<<"Or we can think about setting a lower bound for the strength.."<<endl;
				  return;
				  }
				  */
			}
		}
	};
	void CountingModel::UseAsimovData(int b){  // 0 for background only, 1 for s+b hypothesis
		if(b==0){
			cout<<"		Using Asimov dataset: background only hypothesis"<<endl;
			for(int i=0; i<v_data.size(); i++){
				v_data[i]=0;
				for(int b=v_sigproc[i]; b<vv_exp_sigbkgs.at(i).size(); b++){
					v_data[i]+= vv_exp_sigbkgs.at(i).at(b);
				}
			}
		}else if(b==1){
			cout<<"		Using Asimov dataset: sig+bkg hypothesis"<<endl;
			for(int i=0; i<v_data.size(); i++){
				v_data[i]=0;
				for(int b=0; b<vv_exp_sigbkgs.at(i).size(); b++){
					v_data[i]+= vv_exp_sigbkgs.at(i).at(b);
				}
			}
		}
	}
	CountingModel* CombineModels(CountingModel *cms1, CountingModel *cms2){
		CountingModel *cms = new CountingModel();

		VChannelVSampleVUncertaintyVParameter tmp_vvvv_uncpar = cms1->Get_vvvv_uncpar();
		VChannelVSampleVUncertainty tmp_vvv_pdftype=cms1->Get_vvv_pdftype();	
		VChannelVSampleVUncertainty tmp_vvv_idcorrl=cms1->Get_vvv_idcorrl();	
		vector<string> tmp_v_uncname = cms1->Get_v_uncname();
		for(int ch=0; ch<cms1->NumOfChannels(); ch++){
			//cms.AddChannel(cms1->GetExpectedNumber(ch,0),cms1->GetExpectedNumber(ch,1), cms1->GetExpectedNumber(ch,2), cms1->GetExpectedNumber(ch,3),
			//		cms1->GetExpectedNumber(ch,4), cms1->GetExpectedNumber(ch, 5), cms1->GetExpectedNumber(ch, 6));	
			cms->AddChannel(cms1->GetChannelName(ch), cms1->Get_v_exp_sigbkgs(ch), cms1->GetNSigprocInChannel(ch));
			for(int isamp=0; isamp<tmp_vvv_pdftype.at(ch).size(); isamp++){
				for(int iunc=0; iunc<tmp_vvv_pdftype.at(ch).at(isamp).size(); iunc++){
					if(tmp_vvv_pdftype.at(ch).at(isamp).at(iunc)==typeLogNormal || tmp_vvv_pdftype.at(ch).at(isamp).at(iunc)==typeTruncatedGaussian){
						cms->AddUncertainty(ch, isamp, 
								tmp_vvvv_uncpar.at(ch).at(isamp).at(iunc).at(0), 
								tmp_vvvv_uncpar.at(ch).at(isamp).at(iunc).at(1), 
								tmp_vvv_pdftype.at(ch).at(isamp).at(iunc),
								tmp_v_uncname[tmp_vvv_idcorrl.at(ch).at(isamp).at(iunc)-1]
								);
					}else if(tmp_vvv_pdftype.at(ch).at(isamp).at(iunc)==typeControlSampleInferredLogNormal
							|| tmp_vvv_pdftype.at(ch).at(isamp).at(iunc)==typeGamma){
						cms->AddUncertainty(ch, isamp, 
								tmp_vvvv_uncpar.at(ch).at(isamp).at(iunc).at(0), 
								tmp_vvvv_uncpar.at(ch).at(isamp).at(iunc).at(1), 
								tmp_vvvv_uncpar.at(ch).at(isamp).at(iunc).at(2), 
								tmp_vvv_pdftype.at(ch).at(isamp).at(iunc),
								tmp_v_uncname[tmp_vvv_idcorrl.at(ch).at(isamp).at(iunc)-1]
								);
					}else if(tmp_vvv_pdftype.at(ch).at(isamp).at(iunc)==typeShapeGaussianQuadraticMorph
							|| tmp_vvv_pdftype.at(ch).at(isamp).at(iunc)==typeShapeGaussianLinearMorph){
						int npar = tmp_vvvv_uncpar.at(ch).at(isamp).at(iunc).size();
						cms->AddUncertainty(ch, isamp, 
								npar, &(tmp_vvvv_uncpar.at(ch).at(isamp).at(iunc)[0]), 
								tmp_vvv_pdftype.at(ch).at(isamp).at(iunc),
								tmp_v_uncname[tmp_vvv_idcorrl.at(ch).at(isamp).at(iunc)-1]
								);

					}else {
						cout<<"The pdftype Not implemented yet: "<<tmp_vvv_pdftype.at(ch).at(isamp).at(iunc)<<endl;
					}
				}
			}	
			cms->AddObservedData(ch, cms1->Get_v_data()[ch]);
			cms->SetProcessNames(ch, cms1->GetProcessNames(ch));
		}	

		tmp_vvvv_uncpar = cms2->Get_vvvv_uncpar();
		tmp_vvv_pdftype=cms2->Get_vvv_pdftype();	
		tmp_vvv_idcorrl=cms2->Get_vvv_idcorrl();	
		tmp_v_uncname = cms2->Get_v_uncname();
		for(int ch=0; ch<cms2->NumOfChannels(); ch++){
			int newch = cms->NumOfChannels(); // like ++
			if(cms2->GetDebug()) cout<<"Adding ch = "<<newch<<"th channel from "<<cms2->GetModelName()<<endl;
			cms->AddChannel(cms2->GetChannelName(ch), cms2->Get_v_exp_sigbkgs(ch), cms2->GetNSigprocInChannel(ch));
			if(cms2->GetDebug()) cout<<"now has total "<<cms->NumOfChannels()<<endl;
			for(int isamp=0; isamp<tmp_vvv_pdftype.at(ch).size(); isamp++){
				for(int iunc=0; iunc<tmp_vvv_pdftype.at(ch).at(isamp).size(); iunc++){
					if(tmp_vvv_pdftype.at(ch).at(isamp).at(iunc)==typeLogNormal || tmp_vvv_pdftype.at(ch).at(isamp).at(iunc)==typeTruncatedGaussian){
						cms->AddUncertainty(newch, isamp, 
								tmp_vvvv_uncpar.at(ch).at(isamp).at(iunc).at(0), 
								tmp_vvvv_uncpar.at(ch).at(isamp).at(iunc).at(1), 
								tmp_vvv_pdftype.at(ch).at(isamp).at(iunc),
								tmp_v_uncname[tmp_vvv_idcorrl.at(ch).at(isamp).at(iunc)-1]
								);
					}
					else if(tmp_vvv_pdftype.at(ch).at(isamp).at(iunc)==typeControlSampleInferredLogNormal
							|| tmp_vvv_pdftype.at(ch).at(isamp).at(iunc)==typeGamma){
						cms->AddUncertainty(newch, isamp, 
								tmp_vvvv_uncpar.at(ch).at(isamp).at(iunc).at(0), 
								tmp_vvvv_uncpar.at(ch).at(isamp).at(iunc).at(1), 
								tmp_vvvv_uncpar.at(ch).at(isamp).at(iunc).at(2), 
								tmp_vvv_pdftype.at(ch).at(isamp).at(iunc),
								tmp_v_uncname[tmp_vvv_idcorrl.at(ch).at(isamp).at(iunc)-1]
								);
					}else if(tmp_vvv_pdftype.at(ch).at(isamp).at(iunc)==typeShapeGaussianQuadraticMorph
							|| tmp_vvv_pdftype.at(ch).at(isamp).at(iunc)==typeShapeGaussianLinearMorph){
						int npar = tmp_vvvv_uncpar.at(ch).at(isamp).at(iunc).size();
						cms->AddUncertainty(newch, isamp, 
								npar, &(tmp_vvvv_uncpar.at(ch).at(isamp).at(iunc)[0]), 
								tmp_vvv_pdftype.at(ch).at(isamp).at(iunc),
								tmp_v_uncname[tmp_vvv_idcorrl.at(ch).at(isamp).at(iunc)-1]
								);

					}else {
						cout<<"The pdftype Not implemented yet: "<<tmp_vvv_pdftype.at(ch).at(isamp).at(iunc)<<endl;
					}
				}
			}	
			cms->AddObservedData(newch, cms2->Get_v_data()[ch]);
			cms->SetProcessNames(newch, cms2->GetProcessNames(ch));
		}	

		vector<CountingModel*> ms; ms.push_back(cms1); ms.push_back(cms2);
		for(int m=0; m<ms.size(); m++){
			if(ms[m]->hasParametricShape()){
				RooWorkspace * w = ms[m]->GetWorkSpace();
				vector<int> vsigproc = ms[m]->Get_v_pdfs_sigproc();
				vector< vector<string> > vvpdfs = ms[m]->Get_vv_pdfs();
				vector<TString> vobs = ms[m]->Get_v_pdfs_observables();
				vector< vector<double> > vnorm = ms[m]->Get_vv_pdfs_norm(); 
				vector<string> vchnames = ms[m]->Get_v_pdfs_channelname();
				vector< RooDataSet* > vrds = ms[m]->Get_v_pdfs_roodataset();
				tmp_vvvv_uncpar = ms[m]->Get_vvv_pdfs_normvariation();
				tmp_vvv_pdftype=ms[m]->Get_vvv_pdfs_pdftype();	
				tmp_vvv_idcorrl=ms[m]->Get_vvv_pdfs_idcorrl();	
				tmp_v_uncname = ms[m]->Get_v_uncname();
				vector< vector<double> > vparamunc = ms[m]->Get_v_pdfs_floatParamsUnc(); // from 0 to max_uncorl
				vector<int> vparamIndcorr = ms[m]->Get_v_pdfs_floatParamsIndcorr();      // only for params
				for(int ch=0; ch<vsigproc.size(); ch++){
					int newch = cms->Get_vv_pdfs().size();
					if(ms[m]->GetDebug()) cout<<"Adding ch = "<<newch<<"th channel from "<<ms[m]->GetModelName()<<endl;
					vector<RooAbsPdf*> vs, vb;
					vector<double> vsnorm, vbnorm;
					for(int p =0 ; p<vvpdfs[ch].size(); p++){
						if(p<vsigproc[ch]){
							vs.push_back(w->pdf(vvpdfs[ch][p].c_str()));
							vsnorm.push_back(vnorm[ch][p]);
						}
						else{
							vb.push_back(w->pdf(vvpdfs[ch][p].c_str()));
							vbnorm.push_back(vnorm[ch][p]);
						}
					}
					RooRealVar*x=w->var(vobs[ch]);
					cms->AddChannel(vchnames[ch], x, vs, vsnorm, vb, vbnorm, w);
					if(ms[m]->GetDebug()) cout<<"  AddChannel "<<endl;
					cms->AddObservedDataSet(vchnames[ch], vrds[ch]);
					if(ms[m]->GetDebug()) cout<<"  AddObservedDataSet"<<endl;

					// add uncertainties on norm 
					for(int isamp=0; isamp<vvpdfs[ch].size(); isamp++){
						for(int iunc=0; iunc<tmp_vvv_pdftype.at(ch).at(isamp).size(); iunc++){
							if(tmp_vvv_pdftype.at(ch).at(isamp).at(iunc)==typeLogNormal || tmp_vvv_pdftype.at(ch).at(isamp).at(iunc)==typeTruncatedGaussian){
								cms->AddUncertaintyOnShapeNorm(vchnames[ch], isamp, 
										tmp_vvvv_uncpar.at(ch).at(isamp).at(iunc).at(0), 
										tmp_vvvv_uncpar.at(ch).at(isamp).at(iunc).at(1), 
										tmp_vvv_pdftype.at(ch).at(isamp).at(iunc),
										tmp_v_uncname[tmp_vvv_idcorrl.at(ch).at(isamp).at(iunc)-1]
										);
							}
						}
					}
					if(ms[m]->GetDebug()) cout<<"  AddUncertaintyOnShapeNorm"<<endl;
					if(ms[m]->GetDebug()) cout<<"Added ch = "<<newch<<"th channel from "<<ms[m]->GetModelName()<<endl;
				}

				// add uncertainties on shape parameters
				if(ms[m]->GetDebug()) cout<<" uncname.size="<<tmp_v_uncname.size()<<" vparamunc.size="<<vparamunc.size()<<endl;
				if(ms[m]->GetDebug()) cout<<" "<<vparamIndcorr.size()<<" floating parameters "<<endl;
				for(int ii=0; ii<vparamIndcorr.size(); ii++){
					int id= vparamIndcorr[ii];
					cms->AddUncertaintyOnShapeParam(tmp_v_uncname[id-1], vparamunc[id][0], vparamunc[id][1], vparamunc[id][2], vparamunc[id][3], vparamunc[id][4] );
					if(ms[m]->GetDebug()) cout<<"  AddUncertaintyOnShapeParam "<<tmp_v_uncname[id-1]<<" "
						<<vparamunc[id][0] <<" "
							<<vparamunc[id][1] <<" "
							<<vparamunc[id][2] <<" "
							<<vparamunc[id][3] <<" "
							<<vparamunc[id][4] <<" "
							<<endl;
				}
			}
		}



		return cms;

	}
	double CountingModel::GetExpectedNumber(int index_channel, int index_sample){
		if(index_channel >= vv_exp_sigbkgs_scaled.size()) return -1;
		if(index_sample>= vv_exp_sigbkgs_scaled.at(index_channel).size()) return -1;
		return vv_exp_sigbkgs_scaled.at(index_channel).at(index_sample);
	}
	void CountingModel::RemoveChannelsWithExpectedSignal0orBkg0(int king){
		// should be advoked at the end of model construction, otherwise you will get into trouble ....
		// either before or after "ConfigUncertaintyPdfs"
		std::vector< vector<double> >::iterator iter=vv_exp_sigbkgs.begin(); 
		int skippedchannels_s = 0, skippedchannels_b=0, skippedchannels_sb=0;
		int skippedchannels =0;
		for(; iter!=vv_exp_sigbkgs.end();){
			int position=iter-vv_exp_sigbkgs.begin();
			double bkg = 0, sig=0;
			for(int p=1; p<(*iter).size(); p++){
				bkg += (*iter)[p];
			}
			sig = (*iter)[0];
			if( 
					( (king==1 || king==2) && sig<=0 ) ||
					( (king==0 || king==2) && bkg<=0        ) 
			  ) {
				iter=vv_exp_sigbkgs.erase( iter );
				v_data.erase( v_data.begin()+position );
				v_sigproc.erase( v_sigproc.begin()+position );
				v_channelname.erase(v_channelname.begin()+position);
				vv_exp_sigbkgs_scaled.erase( vv_exp_sigbkgs_scaled.begin()+position );
				vvvv_uncpar.erase( vvvv_uncpar.begin()+position );
				vvv_pdftype.erase( vvv_pdftype.begin()+position );
				vvv_idcorrl.erase( vvv_idcorrl.begin()+position );
				vv_procname.erase( vv_procname.begin()+position );

				if( (king==1||king==2) && sig<=0 )skippedchannels_s++;
				if( (king==0||king==2) && bkg<=0 )skippedchannels_b++;
				if( (king==2) && bkg<=0 && sig<=0 ) skippedchannels_sb++;
				skippedchannels ++;
			}else{
				++iter;
			}	
		}	
		if(king==1 || king==2)	cout<<"\t * "<<skippedchannels_s<<" bins have been removed because of no signal  in those bins *"<<endl;	
		if(king==0 || king==2)	cout<<"\t * "<<skippedchannels_b<<" bins have been removed because of no backgrounds  in those bins *"<<endl;	
		if(king==2)	cout<<"\t * "<<skippedchannels_sb<<" bins have been removed because of no backgrounds and no signal in those bins *"<<endl;	
		cout<<"\t * "<<skippedchannels<<" bins in total were removed *"<<endl;
	}


	// for parametric shapes
	void CountingModel::AddChannel(string channel_name, RooRealVar* observable, vector<RooAbsPdf*> sigPdfs, vector<double> sigNorms, vector<RooAbsPdf*> bkgPdfs, vector<double> bkgNorms, RooWorkspace *w ){
		int signal_processes = sigPdfs.size();
		int bkg_processes = bkgPdfs.size();
		if(signal_processes<=0)  {cout<<"ERROR: you add a channel with number of signal_processes <=0 "<<endl; exit(0);}
		if(bkg_processes<=0)  {cout<<"ERROR: you add a channel with number of bkg_processes <=0 "<<endl; exit(0);}
		v_pdfs_sigproc.push_back(signal_processes);

		if(channel_name==""){
			char tmp[256];
			sprintf(tmp, "channel_%d", v_pdfs_channelname.size());
			channel_name==tmp;
			v_pdfs_channelname.push_back(channel_name);

		}
		else v_pdfs_channelname.push_back(channel_name);


		if(_debug)cout<<" adding channel "<<channel_name<<": "<<sigPdfs.size()<<" signal processes and "<<bkgPdfs.size()<<" bkg processes"<<endl;
		vector<string> vproc; vproc.clear();
		for(int i=0; i<sigPdfs.size(); i++) {
			vproc.push_back(sigPdfs[i]->GetName());
		}
		for(int i=0; i<bkgPdfs.size(); i++) {
			vproc.push_back(bkgPdfs[i]->GetName());
		}
		vv_pdfs_procname.push_back(vproc);

		double tmp_totbkg = 0;

		//RooArgSet* ras = new RooArgSet(*observable);
		//v_pdfs_observables.push_back(ras);

		//observable->SetName(TString(channel_name)+observable->GetName());
		//_w->import(*observable, Rename(TString::Format("%s_%s", channel_name.c_str(), observable->GetName())));
		//_w_varied->import(*observable, Rename(TString::Format("%s_%s", channel_name.c_str(), observable->GetName())));
		_w->import(*observable);
		_w_varied->import(*observable);
		v_pdfs_observables.push_back(observable->GetName());


		vector<string> vsigbkgs; vsigbkgs.clear();
		vector< vector< vector<double> > > vvvuncpar; vvvuncpar.clear();
		vector< vector<int> > vvpdftype; vvpdftype.clear();
		vector< vector<int> > vvidcorrl; vvidcorrl.clear();
		vector<double> vnorms; vnorms.clear();

		vector< vector<double> > vvunc; vvunc.clear();
		vector<int> vpdftype; vpdftype.clear();
		vector<int> vidcorrl; vidcorrl.clear();

		vector<TString> vrrvnorm; vrrvnorm.clear();

		vector<RooRealVar*> vrrvparams; vrrvparams.clear();
		vector< vector<RooRealVar*> > vvrrvparams; vvrrvparams.clear();
		for(int i=0; i<sigPdfs.size(); i++){
			vsigbkgs.push_back(sigPdfs[i]->GetName());
			vvvuncpar.push_back(vvunc);
			vvpdftype.push_back(vpdftype);
			vvidcorrl.push_back(vidcorrl);
			vnorms.push_back(sigNorms[i]);
			TString sn = channel_name; sn+=vproc[i]; sn+="_norm";
			RooRealVar *rrv = new RooRealVar(sn, "", sigNorms[i]);
			vrrvnorm.push_back(sn);
			//rrv->SetName(TString::Format("%s_%s", channel_name.c_str(), rrv->GetName()));
			//sigPdfs[i]->SetName(TString::Format("%s_%s", channel_name.c_str(), sigPdfs[i]->GetName()));
			//_w->import(*rrv, Rename(TString::Format("%s_%s", channel_name.c_str(), rrv->GetName())));
			//_w->import(*sigPdfs[i], Rename(TString::Format("%s_%s", channel_name.c_str(), sigPdfs[i]->GetName())));
			//_w_varied->import(*rrv, Rename(TString::Format("%s_%s", channel_name.c_str(), rrv->GetName())));
			//_w_varied->import(*sigPdfs[i], Rename(TString::Format("%s_%s", channel_name.c_str(), sigPdfs[i]->GetName())));

			_w->import(*rrv);
			_w->import(*sigPdfs[i]);
			_w_varied->import(*rrv);
			_w_varied->import(*sigPdfs[i]);
			RooArgSet *rds	= sigPdfs[i]->getParameters(*observable);
			// need to store the list of parameters and for future modification, fluctuation 
		}

		for(int i=0; i<bkgPdfs.size(); i++){
			vsigbkgs.push_back(bkgPdfs[i]->GetName());
			vvvuncpar.push_back(vvunc);
			vvpdftype.push_back(vpdftype);
			vvidcorrl.push_back(vidcorrl);
			tmp_totbkg+=bkgNorms[i];
			vnorms.push_back(bkgNorms[i]);

			TString sn = channel_name; sn+=vproc[i+sigNorms.size()]; sn+="_norm";
			RooRealVar *rrv = new RooRealVar(sn, "", bkgNorms[i]);
			vrrvnorm.push_back(sn);
			//rrv->SetName(TString::Format("%s_%s", channel_name.c_str(), rrv->GetName()));
			//bkgPdfs[i]->SetName(TString::Format("%s_%s", channel_name.c_str(), bkgPdfs[i]->GetName()));
			//_w->import(*rrv, Rename(TString::Format("%s_%s", channel_name.c_str(), rrv->GetName())));
			//_w->import(*bkgPdfs[i], Rename(TString::Format("%s_%s", channel_name.c_str(), bkgPdfs[i]->GetName())));
			//_w_varied->import(*rrv, Rename(TString::Format("%s_%s", channel_name.c_str(), rrv->GetName())));
			//_w_varied->import(*bkgPdfs[i], Rename(TString::Format("%s_%s", channel_name.c_str(), bkgPdfs[i]->GetName())));
			_w->import(*rrv);
			_w->import(*bkgPdfs[i]);
			_w_varied->import(*rrv);
			_w_varied->import(*bkgPdfs[i]);
		}

		vv_pdfs.push_back(vsigbkgs);
		//vv_exp_sigbkgs_scaled.push_back(vsigbkgs);
		vvv_pdfs_normvariation.push_back(vvvuncpar);
		vvv_pdfs_pdftype.push_back(vvpdftype);
		vvv_pdfs_idcorrl.push_back(vvidcorrl);
		vv_pdfs_normNAME.push_back(vrrvnorm);
		vv_pdfs_norm.push_back(vnorms);
		vv_pdfs_norm_scaled.push_back(vnorms);
		vv_pdfs_norm_varied.push_back(vnorms);
		if(_debug>=100)cout<<"channel "<<channel_name<<": tot.size="<<signal_processes+bkg_processes<<" bkg.size="<<bkg_processes<<endl;


		// construct  model_sb,  model_s,  model_b
		TString s = "SUM::"; s+=channel_name; s+="_sb(";
		for(int i=0; i<vsigbkgs.size(); i++){
			if(i!=0) s+=",";
			s+=vrrvnorm[i]; s+="*"; s+=vsigbkgs[i]; 			
		}
		s+=")";
		_w->factory(s);
		_w_varied->factory(s);
		s = channel_name; s+="_sb";
		v_pdfs_sb.push_back(s);	

		//if(_debug) _w->pdf(s)->Print("V");

		s = "SUM::"; s+=channel_name; s+="_s(";
		for(int i=0; i<signal_processes; i++){
			if(i!=0) s+=",";
			s+=vrrvnorm[i]; s+="*"; s+=vsigbkgs[i]; 			
		}
		s+=")";
		_w->factory(s);
		_w_varied->factory(s);
		s = channel_name; s+="_s";
		v_pdfs_s.push_back(s);	
		//if(_debug) _w->pdf(s)->Print("V");

		s = "SUM::"; s+=channel_name; s+="_b(";
		for(int i=signal_processes; i<vsigbkgs.size(); i++){
			if(i!=signal_processes) s+=",";
			s+=vrrvnorm[i]; s+="*"; s+=vsigbkgs[i]; 			
		}
		s+=")";
		_w->factory(s);
		_w_varied->factory(s);
		s = channel_name; s+="_b";
		v_pdfs_b.push_back(s);	
		//if(_debug) _w->pdf(s)->Print("V");

		RooDataSet * rds = _w->pdf(v_pdfs_b.back()) -> generate(*observable, Extended());
		v_pdfs_roodataset.push_back(rds);
		//v_pdfs_roodataset_toy.push_back(rds);
		//cout<<"H0, number of events generated for channel "<<ch<<": "<<v_pdfs_roodataset_toy[ch]->sumEntries()<<endl;
		if(_debug) 	rds->Print("V");
		vector<double> vdata;
		for(int i=0; i<rds->sumEntries(); i++){
			RooRealVar *r = dynamic_cast<RooRealVar*>(rds->get(i)->first());
			vdata.push_back( r->getVal() );
		}
		vv_pdfs_data.push_back(vdata);


		bHasParametricShape = true;

		if(_debug>=10){
			TString s = "pdf_"; s+= channel_name; s+=".root";
			TFile f(s, "RECREATE") ;
			double xmin = _w->var(v_pdfs_observables.back())->getMin();
			double xmax = _w->var(v_pdfs_observables.back())->getMax();
			TH1F hs("s","s", 100, xmin, xmax);
			TH1F hb("b","b", 100,  xmin, xmax);
			TH1F hsb("sb","sb", 100,  xmin, xmax);
			RooArgSet vars(*(_w->var(v_pdfs_observables.back() ) ) );
			for(int x = 1; x<=100; x++){
				_w->var(v_pdfs_observables.back())->setVal((xmax-xmin)/100.*(x-1)+xmin); 
				hs.SetBinContent(x, _w->pdf(v_pdfs_s.back())->getVal(&vars));
				hb.SetBinContent(x, _w->pdf(v_pdfs_b.back())->getVal(&vars));
				hsb.SetBinContent(x, _w->pdf(v_pdfs_sb.back())->getVal(&vars));
			}

			cout<<"\n model_s"<<endl;
			_w->pdf(v_pdfs_s.back())->getParameters(*(_w->var(v_pdfs_observables.back())))->Print("V");
			cout<<"\n model_b"<<endl;
			_w->pdf(v_pdfs_b.back())->getParameters(*(_w->var(v_pdfs_observables.back())))->Print("V");
			cout<<"\n model_sb"<<endl;
			_w->pdf(v_pdfs_sb.back())->getParameters(*(_w->var(v_pdfs_observables.back())))->Print("V");

			hs.Write();
			hb.Write();
			hsb.Write();
			f.Close();
		}

		//	if(_debug>=10)_w->Print("V");
		if(_debug>=10)_w_varied->Print("V");
	}

	double CountingModel::EvaluateLnQ(int ch, int dataOrToy ){ // 0 for data, 1 for toy
		double ret=0;

		double btot = 0, stot=0;
		for(int i=0; i<vv_pdfs_norm_scaled[ch].size(); i++){
			if(i>=v_pdfs_sigproc[ch]) btot+=vv_pdfs_norm_scaled[ch][i];
			else stot+=vv_pdfs_norm_scaled[ch][i];
		}
		RooArgSet vars(*(_w->var(v_pdfs_observables[ch])));
		if(dataOrToy == 0){
			int ntot = int(v_pdfs_roodataset_toy[ch]->sumEntries());
			for(int i=0; i<v_pdfs_roodataset[ch]->sumEntries(); i++){
				_w->var(v_pdfs_observables[ch])->setVal(( dynamic_cast<RooRealVar*>(v_pdfs_roodataset[ch]->get(i)->first()))->getVal());

				if(_debug>=100){
					if(i==0 or i==ntot-1 or i==ntot/2){
						cout<<"* event "<<i<<":  m= "<<( dynamic_cast<RooRealVar*>(v_pdfs_roodataset[ch]->get(i)->first()))->getVal()<<endl;
						cout<<" pdfs= "<<_w->pdf(v_pdfs_s[ch])->getVal(&vars)<<endl;
						cout<<" pdfb= "<<_w->pdf(v_pdfs_b[ch])->getVal(&vars)<<endl;
						cout<<"stot = "<<stot<<" btot="<<btot<<endl;
					}
				}
				ret+= log(1+ stot/btot*_w->pdf(v_pdfs_s[ch])->getVal(&vars)/_w->pdf(v_pdfs_b[ch])->getVal(&vars));
			}
		}else if(dataOrToy==1){
			int ntot = int(v_pdfs_roodataset_toy[ch]->sumEntries());
			for(int i=0; i<v_pdfs_roodataset_toy[ch]->sumEntries(); i++){
				_w->var(v_pdfs_observables[ch])->setVal(( dynamic_cast<RooRealVar*>(v_pdfs_roodataset_toy[ch]->get(i)->first()))->getVal());
				if(_debug>=100){
					if(i==0 or i==ntot-1 or i==ntot/2){
						cout<<"* event "<<i<<":  m= "<<( dynamic_cast<RooRealVar*>(v_pdfs_roodataset_toy[ch]->get(i)->first()))->getVal()<<endl;
						cout<<" pdfs= "<<_w->pdf(v_pdfs_s[ch])->getVal(&vars)<<endl;
						cout<<" pdfb= "<<_w->pdf(v_pdfs_b[ch])->getVal(&vars)<<endl;
						cout<<"stot = "<<stot<<" btot="<<btot<<endl;
					}
				}
				ret+= log(1+ stot/btot*_w->pdf(v_pdfs_s[ch])->getVal(&vars)/_w->pdf(v_pdfs_b[ch])->getVal(&vars));
			}
		}

		ret-=stot;
		if(_debug>=10)cout<<"EvaluateLnQ of "<<(dataOrToy==0?"data":"toy")<<" in channel ["<<v_pdfs_channelname[ch]<<"]: lnQ= "<<ret<<endl;
		return ret;
	}

	void CountingModel::AddObservedDataSet(int index_channel, RooDataSet* rds){
		//if(v_pdfs_roodataset[index_channel]) delete v_pdfs_roodataset[index_channel];
		int ch = index_channel;

		RooDataSet *tmp = v_pdfs_roodataset[ch];
		v_pdfs_roodataset[ch]=rds;
		delete tmp;

		if(_debug>=10){
			TString s = "data_"; s+= v_pdfs_channelname[ch]; s+=".root";
			TFile f(s, "RECREATE") ;
			TH1F h("data","data", 100, _w->var(v_pdfs_observables[ch])->getMin(), _w->var(v_pdfs_observables[ch])->getMax() );

			for(int i=0; i<v_pdfs_roodataset[ch]->sumEntries(); i++){
				_w->var(v_pdfs_observables[ch])->setVal(( dynamic_cast<RooRealVar*>(v_pdfs_roodataset[ch]->get(i)->first()))->getVal());
				h.Fill(_w->var(v_pdfs_observables[ch])->getVal());
			}
			f.WriteTObject(&h);
			f.Close();
		}
		if(_debug>=10){
			_w->var(v_pdfs_observables[ch])->setVal(6);
			RooArgSet vars(*(_w->var(v_pdfs_observables[ch])));
			cout<<" Norm pdfs= "<<_w->pdf(v_pdfs_s[ch])->getNorm(&vars)<<" val : "<<_w->pdf(v_pdfs_s[ch])->getVal(&vars)<<endl;
			cout<<" Norm pdfb= "<<_w->pdf(v_pdfs_b[ch])->getNorm(&vars)<<" val : "<<_w->pdf(v_pdfs_b[ch])->getVal(&vars)<<endl;
			cout<<" Norm pdfsb= "<<_w->pdf(v_pdfs_sb[ch])->getNorm(&vars)<<" val : "<<_w->pdf(v_pdfs_sb[ch])->getVal(&vars)<<endl;
		}
	}

	void CountingModel::AddObservedDataSet(string channelname, RooDataSet* rds){
		int ch =0;
		for(int c = 0; c<v_pdfs_channelname.size(); c++){
			if(v_pdfs_channelname[c]==channelname) ch=c;
		}
		AddObservedDataSet(ch, rds);
	}
	double CountingModel::EvaluateChi2(double *par){ 
		double ret=0;

		FluctuatedNumbers(par);

		for(int ch=0; ch<vv_pdfs.size(); ch++){
			double btot = 0, stot=0;
			int ntot = int(v_pdfs_roodataset[ch]->sumEntries());
			double tmp=0, retch=0;
			for(int i=0; i<vv_pdfs_norm_varied[ch].size(); i++){
				if(i>=v_pdfs_sigproc[ch]) btot+=vv_pdfs_norm_varied[ch][i];
				else stot+=vv_pdfs_norm_varied[ch][i];
			}
			RooArgSet vars(*(_w_varied->var(v_pdfs_observables[ch])));
			for(int i=0; i<ntot; i++){
				_w_varied->var(v_pdfs_observables[ch])->setVal(( dynamic_cast<RooRealVar*>(v_pdfs_roodataset[ch]->get(i)->first()))->getVal());
				if(_debug>=100){
					if(i==0 or i==ntot-1 or i==ntot/2){
						cout<<"* event "<<i<<":  m= "<<( dynamic_cast<RooRealVar*>(v_pdfs_roodataset[ch]->get(i)->first()))->getVal()<<endl;
						cout<<" pdfs= "<<(par[0]==0?0:_w_varied->pdf(v_pdfs_s[ch])->getVal(&vars))<<endl;
						cout<<" pdfb= "<<_w_varied->pdf(v_pdfs_b[ch])->getVal(&vars)<<endl;
						cout<<" stot = "<<stot<<" btot="<<btot<<endl;
						//tmp = (stot+btot)*_w_varied->pdf(v_pdfs_sb[ch]->GetName())->getVal(&vars);	
						//cout<<" log(event) = "<<log(stot*_w_varied->pdf(v_pdfs_s[ch]->GetName())->getVal(&vars)
						//		+btot*_w_varied->pdf(v_pdfs_b[ch]->GetName())->getVal(&vars))<<endl;
						//cout<<" log(event) ="<<(tmp>0?log(tmp):0)<<endl;
					}
				}
				//tmp = (stot+btot)*_w_varied->pdf(v_pdfs_sb[ch]->GetName())->getVal(&vars);// give some error message ... when r<0
				tmp = 0;  ///////////////
				if(stot!=0) tmp += stot*_w_varied->pdf(v_pdfs_s[ch])->getVal(&vars);  //give some warning message when r=0
				tmp += btot*_w_varied->pdf(v_pdfs_b[ch])->getVal(&vars);

				if(_debug>=100)cout<<" log(event) = "<<log(tmp)<<endl;
				retch -= (tmp>0?log(tmp):0);
			}

			retch+=stot;
			retch+=btot;
			retch-=ntot;
			if(_debug>=10){
				cout<<"EvaluateChi2 in channel ["<<v_pdfs_channelname[ch]<<"]: lnQ= "<<retch<<endl;
				cout<<"\n model_sb"<<endl;
				_w_varied->pdf(v_pdfs_sb[ch])->getParameters(*(_w->var(v_pdfs_observables[ch])))->Print("V");
			}
			ret+=retch;
		}
		return ret;
	}

	double CountingModel::EvaluateGL(int ch, double xr){ 
		double ret=0;

		double btot = 0, stot=0;
		int ntot = int(v_pdfs_roodataset[ch]->sumEntries());
		double tmp;
		for(int i=0; i<vv_pdfs_norm_scaled[ch].size(); i++){
			if(i>=v_pdfs_sigproc[ch]) btot+=vv_pdfs_norm_scaled[ch][i];
			else stot+=vv_pdfs_norm_scaled[ch][i];
		}
		RooArgSet vars(*(_w->var(v_pdfs_observables[ch])));
		for(int i=0; i<ntot; i++){
			_w->var(v_pdfs_observables[ch])->setVal(( dynamic_cast<RooRealVar*>(v_pdfs_roodataset[ch]->get(i)->first()))->getVal());
			if(_debug>=100){
				if(i==0 or i==ntot-1 or i==ntot/2){
					cout<<"* event "<<i<<":  m= "<<( dynamic_cast<RooRealVar*>(v_pdfs_roodataset[ch]->get(i)->first()))->getVal()<<endl;
					cout<<" pdfs= "<<_w->pdf(v_pdfs_s[ch])->getVal(&vars)<<endl;
					cout<<" pdfb= "<<_w->pdf(v_pdfs_b[ch])->getVal(&vars)<<endl;
					cout<<" stot= "<<stot<<" btot="<<btot<<endl;
					cout<<" log(event) = "<<log(stot*_w->pdf(v_pdfs_s[ch])->getVal(&vars)+btot*_w->pdf(v_pdfs_b[ch])->getVal(&vars))<<endl;
				}
			}
			ret+= log( xr*stot*_w->pdf(v_pdfs_s[ch])->getVal(&vars)+btot*_w->pdf(v_pdfs_b[ch])->getVal(&vars));
		}

		if(_debug>=100){
			cout<<"EvaluateGL in channel ["<<v_pdfs_channelname[ch]<<"]: gl = "<<ret<<endl;
			cout<<"\n model_sb"<<endl;
			_w->pdf(v_pdfs_sb[ch])->getParameters(*_w->var(v_pdfs_observables[ch]))->Print("V");
		}
		return ret;
	}
	void CountingModel::SetDataForUnbinned(vector< RooDataSet* > data){
		for(int ch=0; ch<vv_pdfs.size(); ch++){
			delete v_pdfs_roodataset[ch];
		}
		v_pdfs_roodataset.clear();
		for(int ch=0; ch<vv_pdfs.size(); ch++){
			vector<double> vtmp;vtmp.clear();
			v_pdfs_roodataset.push_back(data[ch]);

			for(int i=0; i<v_pdfs_roodataset[ch]->sumEntries(); i++){
				RooRealVar *r = dynamic_cast<RooRealVar*>(v_pdfs_roodataset[ch]->get(i)->first());
				vtmp.push_back( r->getVal() );
			}
			vv_pdfs_data.push_back(vtmp);
		}
	}
	void CountingModel::AddUncertaintyOnShapeNorm(int index_channel, int index_sample, double uncertainty_in_relative_fraction_down, double uncertainty_in_relative_fraction_up, int pdf_type, int index_correlation ){
		// to deal with asymetric uncertainties
		if( uncertainty_in_relative_fraction_down < 0 or uncertainty_in_relative_fraction_up < 0 ) {
			if(pdf_type==typeTruncatedGaussian) {}; //fine
			if( (uncertainty_in_relative_fraction_down <-1 or uncertainty_in_relative_fraction_up <-1) && pdf_type==typeLogNormal) { cout<<"logNormal type uncertainties can't have kappa < 0, exit"<<endl; exit(0);}; //fine
		} 
		if(pdf_type!=typeLogNormal && pdf_type!= typeTruncatedGaussian && pdf_type!=typeGamma ) {
			cout<<"Error: Currently only implemented LogNormal, Gamma and TruncatedGaussian. Your input "<<pdf_type<<" haven't been implemented, exit"<<endl;
			exit(0);
		}
		if(index_correlation <= 0) { 
			cout<<"Error: index_correlation < 0 "<<endl;
			exit(0);
		}
		vector<double> vunc; vunc.clear(); 
		vunc.push_back(uncertainty_in_relative_fraction_down);
		vunc.push_back(uncertainty_in_relative_fraction_up);
		vvv_pdfs_normvariation.at(index_channel).at(index_sample).push_back(vunc);
		vvv_pdfs_pdftype.at(index_channel).at(index_sample).push_back(pdf_type);
		vvv_pdfs_idcorrl.at(index_channel).at(index_sample).push_back(index_correlation);
		//ConfigUncertaintyPdfs();
	}
	void CountingModel::AddUncertaintyOnShapeNorm(int index_channel, int index_sample, double uncertainty_in_relative_fraction_down, double uncertainty_in_relative_fraction_up, int pdf_type, string uncname ){
		int index_correlation = -1; // numeration starts from 1
		for(int i=0; i<v_uncname.size(); i++){
			if(v_uncname[i]==uncname){
				index_correlation = i+1;	
				break;
			}
		}
		if(index_correlation<0)  {
			index_correlation = v_uncname.size()+1;
			v_uncname.push_back(uncname);
		}
		if(_debug>=10) cout<<" AddUncertaintyOnShapeNorm for "<<index_sample<<"th channel, "<<index_sample<<"th proc,  unc ["<<uncname<<"]"<<endl;
		AddUncertaintyOnShapeNorm(index_channel, index_sample, uncertainty_in_relative_fraction_down, uncertainty_in_relative_fraction_up, pdf_type, index_correlation );
	}
	void CountingModel::AddUncertaintyOnShapeNorm(string c, int index_sample, double uncertainty_in_relative_fraction_down, double uncertainty_in_relative_fraction_up, int pdf_type, string uncname ){
		int index_channel = -1;
		for(int i=0; i<v_pdfs_channelname.size(); i++){
			if(v_pdfs_channelname[i]==c) index_channel=i;
		}
		AddUncertaintyOnShapeNorm(index_channel, index_sample, uncertainty_in_relative_fraction_down, uncertainty_in_relative_fraction_up, pdf_type, uncname);

	}

	double CountingModel::EvaluateGL(vector< vector<double> > vnorms, vector<double> vparams, double xr){ 
		double ret=0;
		if(_debug>=10){cout<<"In EvaluateGL:  xr ="<<xr<<endl;};
		for(int i=0; i<v_pdfs_floatParamsName.size(); i++){
			if(_debug>=100) cout<<" EvaluateGL: setting value of parameter ["<<v_pdfs_floatParamsName[i]<<"] = "<<vparams[i]<<endl;
			_w_varied->var(v_pdfs_floatParamsName[i].c_str())->setVal(vparams[i]);
		}
		for(int ch=0; ch<vv_pdfs.size(); ch++){
			double btot = 0, stot=0;
			int ntot = int(v_pdfs_roodataset[ch]->sumEntries());
			double tmp=0;
			for(int i=0; i<vnorms[ch].size(); i++){
				if(i>=v_pdfs_sigproc[ch]) btot+=vnorms[ch][i];
				else stot+=vnorms[ch][i];

				_w_varied->var(vv_pdfs_normNAME[ch][i])->setVal(vnorms[ch][i]);
				if(_debug>=10)cout<<vv_pdfs_normNAME[ch][i]<<" "<<vnorms[ch][i]<<endl;
			}
			RooArgSet vars(*(_w->var(v_pdfs_observables[ch])));
			for(int i=0; i<ntot; i++){
				_w_varied->var(v_pdfs_observables[ch])->setVal(( dynamic_cast<RooRealVar*>(v_pdfs_roodataset[ch]->get(i)->first()))->getVal());
				if(_debug>=100){
					if(i==0 or i==ntot-1 or i==ntot/2){
						cout<<"* event "<<i<<":  m= "<<( dynamic_cast<RooRealVar*>(v_pdfs_roodataset[ch]->get(i)->first()))->getVal()<<endl;
						cout<<" varied_pdfs= "<<_w_varied->pdf(v_pdfs_s[ch])->getVal(&vars)<<endl;
						cout<<" varied_pdfb= "<<_w_varied->pdf(v_pdfs_b[ch])->getVal(&vars)<<endl;
						cout<<" stot= "<<stot<<" btot="<<btot<<endl;
						cout<<" log(event) = "<<log(stot*_w_varied->pdf(v_pdfs_s[ch])->getVal(&vars)
								+btot*_w_varied->pdf(v_pdfs_b[ch])->getVal(&vars))<<endl;
					}
				}
				tmp+= log( xr*stot*_w_varied->pdf(v_pdfs_s[ch])->getVal(&vars)
						+btot*_w_varied->pdf(v_pdfs_b[ch])->getVal(&vars));
			}

			if(_debug>=10){
				cout<<"EvaluateGL in channel ["<<v_pdfs_channelname[ch]<<"]: gl = "<<tmp<<endl;
				cout<<"\n model_sb"<<endl;
				_w_varied->pdf(v_pdfs_sb[ch])->getParameters(*_w->var(v_pdfs_observables[ch]))->Print("V");
			}
			ret+=tmp;
		}
		if(_debug>=10) {
			cout<<" Unbinned Model EvaluateGL all channels = " << ret <<endl;
		}
		return ret;
	}

	void CountingModel::AddUncertaintyOnShapeParam(string pname, double mean, double sigmaL, double sigmaR, double rangeMin, double rangeMax ){
		if(_w_varied->var(pname.c_str())==NULL) {
			cout<<" parameter "<<pname<<" not exist in the added channels,   skip it"<<endl;
			return;
		}

		int index_correlation = -1; // numeration starts from 1
		sigmaL = fabs(sigmaL);
		sigmaR = fabs(sigmaR);
		for(int i=0; i<v_uncname.size(); i++){
			if(v_uncname[i]==pname){
				index_correlation = i+1;	
				break;
			}
		}
		if(index_correlation<0)  {
			index_correlation = v_uncname.size()+1;
			v_uncname.push_back(pname);
		}

		if(rangeMax==rangeMin){
			rangeMin = mean - 4*sigmaL;
			rangeMax = mean + 4*sigmaR;
		}

		vector<double> vunc; vunc.clear(); 
		vunc.push_back(mean);
		vunc.push_back(sigmaL);
		vunc.push_back(sigmaR);
		vunc.push_back(rangeMin);
		vunc.push_back(rangeMax);
		if(v_pdfs_floatParamsUnc.size() <= index_correlation){
			for(int i = v_pdfs_floatParamsUnc.size(); i<=index_correlation; i++){
				v_pdfs_floatParamsUnc.push_back(vunc);
			}
		}else v_pdfs_floatParamsUnc[index_correlation] = vunc;
		if(_debug) cout<<" v_pdfs_floatParamsUnc.size() = "<<v_pdfs_floatParamsUnc.size()<<endl;

		TString s = pname;
		cout<<" * Adding floating parameter: "<<pname<<endl;
		_w_varied->factory(TString::Format("%s_x[%f,%f]", pname.c_str(), rangeMin, rangeMax));
		cout<<pname<<"_x added"<<endl;
		_w_varied->factory(TString::Format("BifurGauss::%s_bfg(%s_x, %f, %f, %f )", pname.c_str(), pname.c_str(), mean, sigmaL, sigmaR));
		cout<<pname<<"_bfg added"<<endl;

		v_pdfs_floatParamsName.push_back(pname);
		v_pdfs_floatParamsIndcorr.push_back(index_correlation);
	}
};
