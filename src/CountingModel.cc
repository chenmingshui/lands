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
//#include <math.h>
using namespace std;
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

	void CountingModel::SetProcessNames(int index_channel, vector<string> vproc){
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
			if(tmpmax>0){
				v_TruncatedGaussian_maxUnc.back()=tmpmax;
			} 
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
	VChannelVSample CountingModel::FluctuatedNumbers(){
		if(_rdm==NULL) {cout<<"Model random gen engine not set yet, exit "<<endl; exit(0);}
		if(!b_systematics) return vv_exp_sigbkgs_scaled;

		double tmp ; 
		vector<double> vrdm; vrdm.clear();
		//if(_debug)cout<<" v_pdftype.size()="<<v_pdftype.size()<<endl;
		for(int i=0; i<v_pdftype.size(); i++){
			//if(_debug)cout<<" vpdftype: "<<i<<"th --> "<<v_pdftype[i]<<endl;
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
					//      one way is to build  TruncatedGaussian function and throw random number from it

					//		vrdm.back()=v_TruncatedGaussian[i]->GetRandom();

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
		return v;
	}
	VIChannel CountingModel::GetToyData_H1(){
		// alternative hypothesis
		VChannelVSample vv = FluctuatedNumbers();

		VIChannel v; v.clear();
		double tmp;
		for(int ch=0; ch<vv.size(); ch++){
			tmp=0;
			for(int isam=0; isam<vv[ch].size(); isam++){
				//start from 0, add signal
				tmp+=vv[ch][isam];	
			}
			v.push_back(_rdm->Poisson(tmp));
			//cout<<"debug tmp="<<tmp<<", pos="<<v[v.size()-1]<<endl;
		}
		//cout<<"delete me, CountingModel::GetToyData_H1 v.size= "<<v.size()<<endl;
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
		cout<<" "<<v_uncname.size()<<" "<<v_pdftype.size()<<endl;
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
	CountingModel CombineModels(CountingModel *cms1, CountingModel *cms2){
		CountingModel cms;

		VChannelVSampleVUncertaintyVParameter tmp_vvvv_uncpar = cms1->Get_vvvv_uncpar();
		VChannelVSampleVUncertainty tmp_vvv_pdftype=cms1->Get_vvv_pdftype();	
		VChannelVSampleVUncertainty tmp_vvv_idcorrl=cms1->Get_vvv_idcorrl();	
		vector<string> tmp_v_uncname = cms1->Get_v_uncname();
		for(int ch=0; ch<cms1->NumOfChannels(); ch++){
			//cms.AddChannel(cms1->GetExpectedNumber(ch,0),cms1->GetExpectedNumber(ch,1), cms1->GetExpectedNumber(ch,2), cms1->GetExpectedNumber(ch,3),
			//		cms1->GetExpectedNumber(ch,4), cms1->GetExpectedNumber(ch, 5), cms1->GetExpectedNumber(ch, 6));	
			cms.AddChannel(cms1->GetChannelName(ch), cms1->Get_v_exp_sigbkgs(ch), cms1->GetNSigprocInChannel(ch));
			for(int isamp=0; isamp<tmp_vvv_pdftype.at(ch).size(); isamp++){
				for(int iunc=0; iunc<tmp_vvv_pdftype.at(ch).at(isamp).size(); iunc++){
					if(tmp_vvv_pdftype.at(ch).at(isamp).at(iunc)==typeLogNormal || tmp_vvv_pdftype.at(ch).at(isamp).at(iunc)==typeTruncatedGaussian){
						cms.AddUncertainty(ch, isamp, 
								tmp_vvvv_uncpar.at(ch).at(isamp).at(iunc).at(0), 
								tmp_vvvv_uncpar.at(ch).at(isamp).at(iunc).at(1), 
								tmp_vvv_pdftype.at(ch).at(isamp).at(iunc),
								tmp_v_uncname[tmp_vvv_idcorrl.at(ch).at(isamp).at(iunc)-1]
								);
					}else if(tmp_vvv_pdftype.at(ch).at(isamp).at(iunc)==typeControlSampleInferredLogNormal
							|| tmp_vvv_pdftype.at(ch).at(isamp).at(iunc)==typeGamma){
						cms.AddUncertainty(ch, isamp, 
								tmp_vvvv_uncpar.at(ch).at(isamp).at(iunc).at(0), 
								tmp_vvvv_uncpar.at(ch).at(isamp).at(iunc).at(1), 
								tmp_vvvv_uncpar.at(ch).at(isamp).at(iunc).at(2), 
								tmp_vvv_pdftype.at(ch).at(isamp).at(iunc),
								tmp_v_uncname[tmp_vvv_idcorrl.at(ch).at(isamp).at(iunc)-1]
								);
					}else if(tmp_vvv_pdftype.at(ch).at(isamp).at(iunc)==typeShapeGaussianQuadraticMorph
							|| tmp_vvv_pdftype.at(ch).at(isamp).at(iunc)==typeShapeGaussianLinearMorph){
						int npar = tmp_vvvv_uncpar.at(ch).at(isamp).at(iunc).size();
						cms.AddUncertainty(ch, isamp, 
								npar, &(tmp_vvvv_uncpar.at(ch).at(isamp).at(iunc)[0]), 
								tmp_vvv_pdftype.at(ch).at(isamp).at(iunc),
								tmp_v_uncname[tmp_vvv_idcorrl.at(ch).at(isamp).at(iunc)-1]
								);

					}else {
						cout<<"The pdftype Not implemented yet: "<<tmp_vvv_pdftype.at(ch).at(isamp).at(iunc)<<endl;
					}
				}
			}	
			cms.AddObservedData(ch, cms1->Get_v_data()[ch]);
			cms.SetProcessNames(ch, cms1->GetProcessNames(ch));
		}	

		tmp_vvvv_uncpar = cms2->Get_vvvv_uncpar();
		tmp_vvv_pdftype=cms2->Get_vvv_pdftype();	
		tmp_vvv_idcorrl=cms2->Get_vvv_idcorrl();	
		tmp_v_uncname = cms2->Get_v_uncname();
		for(int ch=0; ch<cms2->NumOfChannels(); ch++){
			int newch = cms.NumOfChannels(); // like ++
			if(cms2->GetDebug()) cout<<"Adding ch = "<<newch<<"th channel from "<<cms2->GetModelName()<<endl;
			cms.AddChannel(cms2->GetChannelName(ch), cms2->Get_v_exp_sigbkgs(ch), cms2->GetNSigprocInChannel(ch));
			if(cms2->GetDebug()) cout<<"now has total "<<cms.NumOfChannels()<<endl;
			for(int isamp=0; isamp<tmp_vvv_pdftype.at(ch).size(); isamp++){
				for(int iunc=0; iunc<tmp_vvv_pdftype.at(ch).at(isamp).size(); iunc++){
					if(tmp_vvv_pdftype.at(ch).at(isamp).at(iunc)==typeLogNormal || tmp_vvv_pdftype.at(ch).at(isamp).at(iunc)==typeTruncatedGaussian){
						cms.AddUncertainty(newch, isamp, 
								tmp_vvvv_uncpar.at(ch).at(isamp).at(iunc).at(0), 
								tmp_vvvv_uncpar.at(ch).at(isamp).at(iunc).at(1), 
								tmp_vvv_pdftype.at(ch).at(isamp).at(iunc),
								tmp_v_uncname[tmp_vvv_idcorrl.at(ch).at(isamp).at(iunc)-1]
								);
					}
					else if(tmp_vvv_pdftype.at(ch).at(isamp).at(iunc)==typeControlSampleInferredLogNormal
							|| tmp_vvv_pdftype.at(ch).at(isamp).at(iunc)==typeGamma){
						cms.AddUncertainty(newch, isamp, 
								tmp_vvvv_uncpar.at(ch).at(isamp).at(iunc).at(0), 
								tmp_vvvv_uncpar.at(ch).at(isamp).at(iunc).at(1), 
								tmp_vvvv_uncpar.at(ch).at(isamp).at(iunc).at(2), 
								tmp_vvv_pdftype.at(ch).at(isamp).at(iunc),
								tmp_v_uncname[tmp_vvv_idcorrl.at(ch).at(isamp).at(iunc)-1]
								);
					}else if(tmp_vvv_pdftype.at(ch).at(isamp).at(iunc)==typeShapeGaussianQuadraticMorph
							|| tmp_vvv_pdftype.at(ch).at(isamp).at(iunc)==typeShapeGaussianLinearMorph){
						int npar = tmp_vvvv_uncpar.at(ch).at(isamp).at(iunc).size();
						cms.AddUncertainty(newch, isamp, 
								npar, &(tmp_vvvv_uncpar.at(ch).at(isamp).at(iunc)[0]), 
								tmp_vvv_pdftype.at(ch).at(isamp).at(iunc),
								tmp_v_uncname[tmp_vvv_idcorrl.at(ch).at(isamp).at(iunc)-1]
								);

					}else {
						cout<<"The pdftype Not implemented yet: "<<tmp_vvv_pdftype.at(ch).at(isamp).at(iunc)<<endl;
					}
				}
			}	
			cms.AddObservedData(newch, cms2->Get_v_data()[ch]);
			cms.SetProcessNames(newch, cms2->GetProcessNames(ch));
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
};
