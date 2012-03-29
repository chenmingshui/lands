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
#include "UtilsROOT.h"
#include <iostream>
#include <cstdio>
#include <cmath>
#include <memory>
#include "TF1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h"
#include "RooAbsPdf.h"
#include "RooAbsData.h"
#include "RooDataHist.h"
#include "RooDataSet.h"
#include "RooRealVar.h"
#include "RooBifurGauss.h"
#include "RooWorkspace.h"
#include "RooHistFunc.h"
#include "CLsLimit.h"
#include "TGraphAsymmErrors.h"
//#include "RooStats/RooStatsUtils.h"

using namespace std;
using namespace RooFit;
//using namespace RooStats;
namespace lands{
    CountingModel::CountingModel(){
        v_data.clear();
        v_data_real.clear();
        vv_exp_sigbkgs.clear();
        vv_exp_sigbkgs_scaled.clear();
        vv_randomized_sigbkgs.clear();
        vv_randomized_sigbkgs_scaled.clear();
        vv_fitted_sigbkgs.clear();
        vv_fitted_sigbkgs_scaled.clear();
        vvvv_uncpar.clear();
        vvv_pdftype.clear();
        vvv_idcorrl.clear();
        v_channelname.clear();
        _rdm=0;
        b_systematics=0;
        b_ForceSymmetryError=0;
        b_MultiSigProcShareSamePDF=0;
        v_TruncatedGaussian_maxUnc.clear();
        v_pdftype.clear();
        _common_signal_strength=1;
        max_uncorrelation=0;
        b_AllowNegativeSignalStrength = 1;
        v_GammaN.clear();
        v_uncname.clear();
        v_uncFloatInFit.clear();
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
        vv_pdfs_norm_randomized.clear();	
        vv_pdfs_norm_randomized_scaled.clear();	
        vv_pdfs_norm_fitted.clear();	
        vv_pdfs_norm_fitted_scaled.clear();	
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
        v_pdfs_roodataset_real.clear();
        v_pdfs_roodataset_tmp.clear();
        v_pdfs_roodataset_asimovb.clear();
        v_pdfs_roodataset_asimovsb.clear();
        vv_pdfs_normNAME.clear();
        vv_pdfs_extranormNAME.clear();
        vvv_pdfs_paramsRRV.clear();

        v_pdfs_floatParamsName.clear();
        v_pdfs_floatParamsIndcorr.clear();
        v_pdfs_floatParamsType.clear();
        v_pdfs_floatParamsUnc.clear();


        _w = new RooWorkspace("w");
        _w_varied = new RooWorkspace();


        _fittedParsInData_bonly = 0;
        _fittedParsInData_sb= 0;
        _fittedParsInData_global= 0;
        _fittedParsInPseudoData_bonly = 0;
        _fittedParsInPseudoData_sb= 0;
        _fittedParsInPseudoData_global= 0;
        _tossToyConvention = 0;
        _UseBestEstimateToCalcQ = 1;

        map_param_sources.clear();

        _norminalPars = 0;
        _norminalParsSigma = 0;
        _randomizedPars= 0;

        vv_pdfs_statusUpdated.clear();
        vv_statusUpdated.clear();
        vvp_connectNuisBinProc.clear();
        vvp_pdfs_connectNuisBinProc.clear();

	_bFixingPOIs = false;
	vsPOIsToBeFixed.clear();


        minuitSTRATEGY = 0; 
        maximumFunctionCallsInAFit=5000;
	minuitTolerance = 0.001;
	nuisancesRange = 5.;

        _PhysicsModel = typeStandardModel;
    }
    CountingModel::~CountingModel(){
        v_data.clear();
        v_data_real.clear();
        vv_exp_sigbkgs.clear();
        vv_exp_sigbkgs_scaled.clear();
        vv_randomized_sigbkgs.clear();
        vv_randomized_sigbkgs_scaled.clear();
        vv_fitted_sigbkgs.clear();
        vv_fitted_sigbkgs_scaled.clear();
        vvvv_uncpar.clear();
        vvv_pdftype.clear();
        vvv_idcorrl.clear();
        v_channelname.clear();
        _rdm=0;

        v_TruncatedGaussian_maxUnc.clear();
        v_pdftype.clear();
        v_GammaN.clear();
        v_uncname.clear();
        v_uncFloatInFit.clear();
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
        vv_pdfs_norm_randomized.clear();	
        vv_pdfs_norm_randomized_scaled.clear();	
        vv_pdfs_norm_fitted.clear();	
        vv_pdfs_norm_fitted_scaled.clear();	
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
        v_pdfs_roodataset_real.clear();
        v_pdfs_roodataset_tmp.clear();
        v_pdfs_roodataset_asimovb.clear();
        v_pdfs_roodataset_asimovsb.clear();
        vv_pdfs_normNAME.clear();
        vv_pdfs_extranormNAME.clear();
        vvv_pdfs_paramsRRV.clear();
        v_pdfs_floatParamsName.clear();
        v_pdfs_floatParamsIndcorr.clear();
        v_pdfs_floatParamsType.clear();
        v_pdfs_floatParamsUnc.clear();

        /*
           if(_fittedParsInData_bonly) delete [] _fittedParsInData_bonly;
           if(_fittedParsInData_sb) delete [] _fittedParsInData_sb;
           if(_fittedParsInData_global) delete [] _fittedParsInData_global;
           if(_fittedParsInPseudoData_bonly) delete [] _fittedParsInPseudoData_bonly;
           if(_fittedParsInPseudoData_sb) delete [] _fittedParsInPseudoData_sb;
           if(_fittedParsInPseudoData_global) delete [] _fittedParsInPseudoData_global;
           */

        _fittedParsInData_bonly=0;
        _fittedParsInData_sb=0;
        _fittedParsInData_global=0;
        _fittedParsInPseudoData_bonly=0;
        _fittedParsInPseudoData_sb=0;
        _fittedParsInPseudoData_global=0;

        map_param_sources.clear();

        vv_pdfs_statusUpdated.clear();
        vv_statusUpdated.clear();
        vvp_connectNuisBinProc.clear();
        vvp_pdfs_connectNuisBinProc.clear();

        delete _w;
        delete _w_varied;

        if(_norminalPars) delete _norminalPars;
        if(_norminalParsSigma) delete _norminalParsSigma;
        if(_randomizedPars) delete _randomizedPars;
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
        v_data_real.push_back(tmp_totbkg);	
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
        v_data_real.push_back(tmp_totbkg);	
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
    void CountingModel::TagUncertaintyFloatInFit(string uncname, bool b){
        bool tagged = false;
        for(int i=0; i<v_uncname.size(); i++){
            if(v_uncname[i]==uncname){
                if(v_uncFloatInFit[i]==false and b==true) { cout<<"ERROR: ["<<uncname<<"]  has been tagged [nofloat] before, and there is conflict now"<<endl; exit(1); }
                v_uncFloatInFit[i]=b;
                tagged = true;
            }
        }	
        if(!tagged) {cout<<"ERROR: in TagUncertaintyFloatInFit, "<<uncname<<" not found, may not added yet "<<endl;
            for(int i=0; i<v_uncname.size(); i++) cout<<v_uncname[i]<<endl;
            exit(1);}
    }
    void CountingModel::TagUncertaintyFloatInFit(int i, bool b){
        if(i>=v_uncFloatInFit.size() or i<0) { cout<<"ERROR TagUncertaintyFloatInFit index "<<i<<" out of range "<<endl; exit(1); }
        if(v_uncFloatInFit[i]==false and b==true) { cout<<"ERROR: ["<<i<<"]  has been tagged [nofloat] before, and there is conflict now"<<endl; exit(1); }
        v_uncFloatInFit[i]=b;	
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
            v_uncFloatInFit.push_back(1);
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
            v_uncFloatInFit.push_back(1);
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
            v_uncFloatInFit.push_back(1);
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
        if(pdf_type!=typeLogNormal && pdf_type!= typeTruncatedGaussian && pdf_type!=typeGamma && pdf_type!=typeFlat) {
            cout<<"Error: Currently only implemented LogNormal, Gamma and TruncatedGaussian, flat. Your input "<<pdf_type<<" haven't been implemented, exit"<<endl;
            exit(0);
        }
        if(index_correlation <= 0) { 
            cout<<"Error: index_correlation < 0 "<<endl;
            exit(0);
        }
        vector<double> vunc; vunc.clear(); 
        if(b_ForceSymmetryError) {
            if (uncertainty_in_relative_fraction_up > uncertainty_in_relative_fraction_down) uncertainty_in_relative_fraction_down=uncertainty_in_relative_fraction_up;
            if(uncertainty_in_relative_fraction_down > uncertainty_in_relative_fraction_up ) uncertainty_in_relative_fraction_up = uncertainty_in_relative_fraction_down;
        }
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
        v_data_real.at(index_channel)=num_data;
    }
    void CountingModel::AddObservedData(string c, double num_data){
        int index_channel = -1;
        for(int i=0; i<v_channelname.size(); i++){
            if(v_channelname[i]==c) index_channel=i;
        }
        v_data.at(index_channel)=num_data;
        v_data_real.at(index_channel)=num_data;
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

        vvp_connectNuisBinProc.clear();
        vvp_pdfs_connectNuisBinProc.clear();

        for(int ch=0; ch<vvv_idcorrl.size(); ch++){
            for(int isam=0; isam<vvv_idcorrl.at(ch).size(); isam++){
                for(int iunc=0; iunc<vvv_idcorrl.at(ch).at(isam).size(); iunc++){
                    if(_debug>=100)	 cout<< "in counting ConfigUncertaintyPdfs: ch "<<ch<<" isam "<<isam<<" iunc "<<iunc<<endl;
                    int indexcorrl = vvv_idcorrl.at(ch).at(isam).at(iunc);
                    if(max_uncorrelation<indexcorrl) max_uncorrelation=indexcorrl;
                }
            }
        }
        for(int ch=0; ch<vvv_pdfs_idcorrl.size(); ch++){
            for(int isam=0; isam<vvv_pdfs_idcorrl.at(ch).size(); isam++){
                for(int iunc=0; iunc<vvv_pdfs_idcorrl.at(ch).at(isam).size(); iunc++){
                    if(_debug>=100)	 cout<< "in shape ConfigUncertaintyPdfs: ch "<<ch<<" isam "<<isam<<" iunc "<<iunc<<endl;
                    int indexcorrl = vvv_pdfs_idcorrl.at(ch).at(isam).at(iunc);
                    if(max_uncorrelation<indexcorrl) max_uncorrelation=indexcorrl;
                }
            }
        }
        for(int i=0; i<v_pdfs_floatParamsIndcorr.size(); i++){
            if(max_uncorrelation < v_pdfs_floatParamsIndcorr[i]) max_uncorrelation=v_pdfs_floatParamsIndcorr[i];
        }

        if(max_uncorrelation < v_uncname.size()) max_uncorrelation = v_uncname.size();

        vvp_connectNuisBinProc.clear();
        vvp_pdfs_connectNuisBinProc.clear();
        vvp_connectNuisBinProc.resize(max_uncorrelation+1);
        vvp_pdfs_connectNuisBinProc.resize(max_uncorrelation+1);// including the signal strength

        for(int i=0; i<=max_uncorrelation; i++){
            v_TruncatedGaussian_maxUnc.push_back(-1);
            v_pdftype.push_back(-1);	
            v_GammaN.push_back(-1);
            double tmpmax=-1;
            for(int ch=0; ch<vvv_idcorrl.size(); ch++){
                for(int isam=0; isam<vvv_idcorrl.at(ch).size(); isam++){
                    if(i==0 && isam<v_sigproc[ch])vvp_connectNuisBinProc[0].push_back(std::make_pair(ch, isam));
                    for(int iunc=0; iunc<vvv_idcorrl.at(ch).at(isam).size(); iunc++){
//cout<<" *** DELETEME ConfigUncertaintyPdfs,  vvv_pdftype[ch][isam][iunc]="<<vvv_pdftype.at(ch).at(isam).at(iunc)<<endl;
                        if(_debug>=101)	 cout<< "in counting ConfigUncertaintyPdfs: ch "<<ch<<" isam "<<isam<<" iunc "<<iunc<<endl;
                        int indexcorrl = vvv_idcorrl.at(ch).at(isam).at(iunc);
                        if(indexcorrl==i && v_pdftype.back()<0 ){
                            v_pdftype.back()=vvv_pdftype.at(ch).at(isam).at(iunc);
//			   cout<<" *** DELETEME ConfigUncertaintyPdfs, v_pdftype.back = "<<v_pdftype.back()<<endl;;
                        }

                        if(indexcorrl==i)vvp_connectNuisBinProc[indexcorrl].push_back(std::make_pair(ch, isam));

                        if(indexcorrl==i && vvv_pdftype.at(ch).at(isam).at(iunc)== typeTruncatedGaussian ){
                            if(tmpmax< fabs(vvvv_uncpar.at(ch).at(isam).at(iunc).at(0)) ) tmpmax=fabs(vvvv_uncpar.at(ch).at(isam).at(iunc).at(0));	
                            if(tmpmax< fabs(vvvv_uncpar.at(ch).at(isam).at(iunc).at(1)) ) tmpmax=fabs(vvvv_uncpar.at(ch).at(isam).at(iunc).at(1));	
                        } 
                        if(indexcorrl==i && vvv_pdftype.at(ch).at(isam).at(iunc)== typeGamma){
                            if(isam<v_sigproc[ch])vvp_connectNuisBinProc[0].push_back(std::make_pair(ch,isam));
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
                            if( v_pdftype.back()!=vvv_pdftype.at(ch).at(isam).at(iunc)
                                    && ( vvv_pdftype.at(ch).at(isam).at(iunc)!=typeShapeGaussianQuadraticMorph &&
                                        vvv_pdftype.at(ch).at(isam).at(iunc)!=typeShapeGaussianLinearMorph) 
                                    && v_pdftype.back()!=typeShapeGaussianLinearMorph 
                                    && v_pdftype.back()!=typeShapeGaussianQuadraticMorph
                              ){
                                cout<<" Error:  two uncertainties with 100% correlation must be with same pdftype, exit "<<endl;
                                cout<<" Independant unc "<<indexcorrl<<"th, name "<<v_uncname[indexcorrl-1];
                                cout<<" should be "<<v_pdftype.back()<<", but "<<vvv_pdftype.at(ch).at(isam).at(iunc)<<" for ch"<<ch<<"("<<v_channelname[ch]<<")"
                                    <<" isam"<<isam<<"("<<vv_procname[ch][isam]<<")"
                                    <<" iunc"<<iunc<<"("<<v_uncname[indexcorrl-1]<<")"<<endl;
                                cout<<" The conflict unc at above process is for "<<vvv_idcorrl[ch][isam][iunc]
                                    <<"th source, name "<<v_uncname[vvv_idcorrl[ch][isam][iunc]-1]<<endl;
                                exit(0);
                            }
                        }
                    }
                }
            }
            for(int ch=0; ch<vvv_pdfs_idcorrl.size(); ch++){
                if(_debug>=100)cout<<ch<<endl;
                for(int isam=0; isam<vvv_pdfs_idcorrl.at(ch).size(); isam++){
                    if(_debug>=100)cout<<isam<<endl;
                    for(int iunc=0; iunc<vvv_pdfs_idcorrl.at(ch).at(isam).size(); iunc++){
                        if(_debug>=100)cout<<iunc<<endl;
                        if(_debug>=100)	 cout<< "in shape ConfigUncertaintyPdfs: ch "<<v_pdfs_channelname[ch]<<" isam "<<vv_pdfs[ch][isam]<<endl;
                        int indexcorrl = vvv_pdfs_idcorrl.at(ch).at(isam).at(iunc);
                        if(_debug>=10) cout<<"      "<<v_uncname[indexcorrl-1]<<endl;
                        if(indexcorrl==i && v_pdftype.back()<0 ){
                            v_pdftype.back()=vvv_pdfs_pdftype.at(ch).at(isam).at(iunc);
                        }
                        if(indexcorrl==i && vvv_pdfs_pdftype.at(ch).at(isam).at(iunc)== typeTruncatedGaussian ){
                            if(tmpmax< fabs(vvv_pdfs_normvariation.at(ch).at(isam).at(iunc).at(0)) ) tmpmax=fabs(vvv_pdfs_normvariation.at(ch).at(isam).at(iunc).at(0));	
                            if(tmpmax< fabs(vvv_pdfs_normvariation.at(ch).at(isam).at(iunc).at(1)) ) tmpmax=fabs(vvv_pdfs_normvariation.at(ch).at(isam).at(iunc).at(1));	
                        } 
                        if(indexcorrl==i && vvv_pdfs_pdftype.at(ch).at(isam).at(iunc)== typeGamma){
                            if(vvv_pdfs_normvariation.at(ch).at(isam).at(iunc).at(0)>0 && 
                                    fabs(vvv_pdfs_normvariation.at(ch).at(isam).at(iunc).at(0)*vvv_pdfs_normvariation.at(ch).at(isam).at(iunc).at(2) - vv_pdfs_norm_varied.at(ch).at(isam)) / vvv_pdfs_normvariation.at(ch).at(isam).at(iunc).at(2)/vvv_pdfs_normvariation.at(ch).at(isam).at(iunc).at(0)>0.2
                              ) {
                                // if there is extra norm inside workspace, then probably will not have gmN uncertainty 
                                cout<<"SHOULD BE FATAL: Shape channel "<<ch<<"th, "<<v_pdfs_channelname[ch]<<": process "<<isam<<" using gamma pdf, but rho*B!=b, please check"<<endl; 
                                cout<<"rho="<<vvv_pdfs_normvariation[ch][isam][iunc][0]<<"  B="<<vvv_pdfs_normvariation[ch][isam][iunc][2]<<" b="<<vv_pdfs_norm[ch][isam]<<endl;
                                cout<<"to proceed, rho are modified. "<<endl;
                                vvv_pdfs_normvariation[ch][isam][iunc][0] = vv_pdfs_norm[ch][isam]/vvv_pdfs_normvariation[ch][isam][iunc][2];
                                //exit(0);
                            }	
                            v_GammaN.back()=vvv_pdfs_normvariation.at(ch).at(isam).at(iunc).at(2)+1;	

                        }
                        if(indexcorrl==i && v_pdftype.back()>0 ){
                            if( v_pdftype.back()!=vvv_pdfs_pdftype.at(ch).at(isam).at(iunc)
                                    && ( vvv_pdfs_pdftype.at(ch).at(isam).at(iunc)!=typeShapeGaussianQuadraticMorph &&
                                        vvv_pdfs_pdftype.at(ch).at(isam).at(iunc)!=typeShapeGaussianLinearMorph) 
                                    && v_pdftype.back()!=typeShapeGaussianLinearMorph 
                                    && v_pdftype.back()!=typeShapeGaussianQuadraticMorph
                              ){
                                cout<<" Error:  two uncertainties with 100% correlation must be with same pdftype, exit "<<endl;
                                cout<<" Independant unc "<<indexcorrl<<"th, name "<<v_uncname[indexcorrl-1];
                                cout<<" should be "<<v_pdftype.back()<<", but "<<vvv_pdfs_pdftype.at(ch).at(isam).at(iunc)<<" for shape ch"<<ch<<"("<<v_pdfs_channelname[ch]<<")"
                                    <<" isam"<<isam<<"("<<vv_pdfs[ch][isam]<<")"
                                    <<" iunc"<<iunc<<"("<<v_uncname[indexcorrl-1]<<")"<<endl;
                                cout<<" The conflict unc at above process is for "<<vvv_pdfs_idcorrl[ch][isam][iunc]
                                    <<"th source, name "<<v_uncname[vvv_pdfs_idcorrl[ch][isam][iunc]-1]<<endl;
                                exit(0);
                            }
                        }
                        if(_debug>=100) cout<<ch<<" "<<isam<< " "<<iunc<<endl;
                    }
                }
            }
            if(tmpmax>0){
                v_TruncatedGaussian_maxUnc.back()=tmpmax;
            } 
        }

        vector< vector<RooArgSet* > > vvr;
        //vvr.resize(vv_pdfs.size());
        for(int ch = 0; ch<vv_pdfs.size(); ch++ ){
            vector< RooArgSet* > vr; vr.clear();
            //vvr[ch].resize(vv_pdfs[ch].size());
            for(int p=0; p<vv_pdfs[ch].size(); p++){
                //RooArgSet* params = _w_varied->pdf(vv_pdfs[ch][p].c_str())->getParameters(*(_w->var(v_pdfs_observables[ch])));
                RooArgSet* params = _w_varied->pdf(vv_pdfs[ch][p].c_str())->getParameters(*(v_pdfs_roodataset_real[ch]->get()));
                //vvr[ch][p]=params;
                vr.push_back(params);

                //	cout<<"DELETEME pdf: "<<vv_pdfs[ch][p]<<" containes: "<<endl;
                //	params->Print("V");
            }
            vvr.push_back(vr);
        }


	vvv_pdfs_nuisancesindex.clear();
	for(int ch = 0; ch<vv_pdfs.size(); ch++ ){
		vector< vector<int> > vv;
		for(int p=0; p<vv_pdfs[ch].size(); p++){
			vector< int > v;
			vv.push_back(v);
		}
		vvv_pdfs_nuisancesindex.push_back(vv);
	}
	for(int i=0; i<v_pdfs_floatParamsName.size(); i++){
            //v_pdftype[v_pdfs_floatParamsIndcorr[i]] = typeBifurcatedGaussian;
            v_pdftype[v_pdfs_floatParamsIndcorr[i]] = v_pdfs_floatParamsType[i];
            // identify here which pdfs are affected by this parameter
            // and update the vvp_pdfs_connectNuisBinProc
	
	int nuiInd = -1;
        for(int j=0; j<v_uncname.size(); j++){
		if(v_uncname[j]==v_pdfs_floatParamsName[i]) nuiInd=j+1;
	}
	if(nuiInd<0){ cout<<" "<<v_pdfs_floatParamsName[i]<<" doesn't appear in uncname list" <<endl; exit(1); }
            //cout<<" DELETEME floatparam "<<v_pdfs_floatParamsName[i]<<endl;
            for(int ch = 0; ch<vv_pdfs.size(); ch++ ){
                for(int p=0; p<vv_pdfs[ch].size(); p++){
                    std::auto_ptr<TIterator> iter(vvr[ch][p]->createIterator());
                    //		cout<<" pdf: "<<vv_pdfs[ch][p]<<endl;
                    for (RooAbsArg *par = (RooAbsArg *) iter->Next(); par != 0; par = (RooAbsArg *) iter->Next()) {
                        //			cout<<" DELETEME par "<<par->GetName()<<endl;
                        if(par->GetName()==v_pdfs_floatParamsName[i]) {
				vvp_pdfs_connectNuisBinProc[v_pdfs_floatParamsIndcorr[i]].push_back(make_pair(ch,p));
				vvv_pdfs_nuisancesindex[ch][p].push_back(nuiInd);
				continue;
			}
                    }
                }
            }
        }
        for(int i=0; i<v_uncname.size(); i++){
            if(v_pdftype[i+1]==-1) v_pdftype[i+1]=typeLogNormal; // some uncertainty source not affect normalization but shape 
            // put here, just because need it for dummy purpose, --> generate unit gaussian random number 
            //  FIXME 
            // need to be able to identify which parameters are affected by this uncertainty source here 
            // and update the vvp_pdfs_connectNuisBinProc

	    //cout<<"DELETEME v_pdftype[i+1]="<<v_pdftype[i+1]<<endl;
        }
        MakeListOfShapeUncertainties();

        _norminalPars = new double[max_uncorrelation+1];
        _norminalParsSigma = new double[max_uncorrelation+1];
        _randomizedPars= new double[max_uncorrelation+1];
	_norminalPars[0]=0;
	_norminalParsSigma[0]=-9999;
        for(int i=1; i<=max_uncorrelation; i++){
            if(_debug)cout<<" ** DELETEME wrong"<<i<<" "<<v_uncname[i-1]<<" v_pdftype = "<<v_pdftype[i]<<endl;
            switch (v_pdftype[i]){
                case typeShapeGaussianLinearMorph:
                case typeShapeGaussianQuadraticMorph:
                case typeLogNormal:
                case typeTruncatedGaussian :
                    _norminalPars[i]=0;
                    _norminalParsSigma[i]=1;
                    break;
                case typeGamma:
                    _norminalPars[i]=v_GammaN[i];
                    _norminalParsSigma[i]=sqrt(v_GammaN[i]);
                    break;
                case typeBifurcatedGaussian:
                    _norminalPars[i]=v_pdfs_floatParamsUnc[i][0];
                    _norminalParsSigma[i]= (v_pdfs_floatParamsUnc[i][1]>v_pdfs_floatParamsUnc[i][2]?v_pdfs_floatParamsUnc[i][1]:v_pdfs_floatParamsUnc[i][2]);
                    break;
                case typeFlat:
                    if(_w->var(v_uncname[i-1].c_str())!=NULL){
                        _norminalPars[i]=(v_pdfs_floatParamsUnc[i][0]-v_pdfs_floatParamsUnc[i][3])/v_pdfs_floatParamsUnc[i][1];
                        if(_debug)cout<<" DELETEME flatParam "<<": norminal="<<v_pdfs_floatParamsUnc[i][0]
                            <<" max-min="<<v_pdfs_floatParamsUnc[i][1]<<" min="
                                <<v_pdfs_floatParamsUnc[i][3]<<" max="<<v_pdfs_floatParamsUnc[i][4]<<endl;
                    }else _norminalPars[i] = 0.5;
		    _norminalParsSigma[i] = -9999.; // it's flat ...

                    break;
                default:
                    cout<<"pdftype not yet defined:  "<<v_pdftype[i]<<endl;
                    cout<<"**********"<<endl;
                    exit(0);
            }
            if(_debug)	 cout<<" DELETEME _norminalPars "<<_norminalPars[i]<<endl; 
        }

        for(int i=0; i<vvp_pdfs_connectNuisBinProc.size(); i++){
            if(_debug>=10)cout<<"parameter "<<(i==0?"r":v_uncname[i-1])<<" affects following ch/p "<<endl;
            for(int j=0; j<vvp_pdfs_connectNuisBinProc[i].size(); j++)
                if(_debug)cout<<" channel: "<<v_pdfs_channelname[vvp_pdfs_connectNuisBinProc[i][j].first]<<" process: "<<vv_pdfs[vvp_pdfs_connectNuisBinProc[i][j].first][vvp_pdfs_connectNuisBinProc[i][j].second]<<endl;;
            for(int j=0; j<vvp_connectNuisBinProc[i].size(); j++)
                if(_debug>=10)	cout<<" channel: "<<v_channelname[vvp_connectNuisBinProc[i][j].first]<<" process: "<<vv_procname[vvp_connectNuisBinProc[i][j].first][vvp_connectNuisBinProc[i][j].second]<<endl;;
        }

    }	
    void CountingModel::MakeListOfShapeUncertainties(){
        vvv_shapeuncindex.clear();
        for(int ch=0; ch<vvv_idcorrl.size(); ch++){
            vector< vector<int> > vvp; vvp.clear();
            for(int isam=0; isam<vvv_idcorrl.at(ch).size(); isam++){
                vector<int> vshape; vshape.clear();
                for(int iunc=0; iunc<vvv_idcorrl.at(ch).at(isam).size(); iunc++){
                    int indexcorrl = vvv_idcorrl.at(ch).at(isam).at(iunc);
                    if(v_pdftype[indexcorrl]==typeShapeGaussianQuadraticMorph or v_pdftype[indexcorrl]==typeShapeGaussianLinearMorph){
                        vshape.push_back(iunc);
                    }
                }
                vvp.push_back(vshape);
            }
            vvv_shapeuncindex.push_back(vvp);
        }

    }
    VChannelVSample CountingModel::FluctuatedNumbers(double *par, bool scaled, int bUseBestEstimateToCalcQ, bool includeCountingParts){
        RooRealVar *tmprrv;
        // FIXME  need to think about what's vv_pdfs_norm_varied,  where it's changed,  and whether to replace it with vv_pdfs_norm_randomized_scaled .... ? 
        if(_rdm==NULL) {cout<<"Model random gen engine not set yet, exit "<<endl; exit(0);}
        if(!b_systematics) {
            if(bHasParametricShape){
                //vv_pdfs_params_varied = bUseBestEstimateToCalcQ?vv_pdfs_params:vv_pdfs_params_randomized;
                if(bUseBestEstimateToCalcQ==1)vv_pdfs_norm_varied = scaled?vv_pdfs_norm_scaled:vv_pdfs_norm;  //FIXME HGG need account for extraNorm
                else if(bUseBestEstimateToCalcQ==0) vv_pdfs_norm_varied = scaled?vv_pdfs_norm_randomized_scaled:vv_pdfs_norm_randomized;
                else vv_pdfs_norm_varied = scaled?vv_pdfs_norm_fitted_scaled:vv_pdfs_norm_fitted;
                for(int ch=0; ch<vv_pdfs.size(); ch++){
                    int nsigproc = v_pdfs_sigproc[ch];
                    for(int isam=0; isam<vv_pdfs[ch].size(); isam++){
                        _w_varied->var(vv_pdfs_normNAME[ch][isam])->setVal(vv_pdfs_norm_varied[ch][isam]);
                    }
                }
            }
            if(bUseBestEstimateToCalcQ==1)return scaled?vv_exp_sigbkgs_scaled:vv_exp_sigbkgs;
            else if(bUseBestEstimateToCalcQ==0) return scaled?vv_randomized_sigbkgs_scaled:vv_randomized_sigbkgs;
            else  return scaled?vv_fitted_sigbkgs_scaled:vv_fitted_sigbkgs; // scaled?fittedNusancesInSB:fittedNusancesInBonly
        }

        double tmp ; 
        vector<double> vrdm; vrdm.clear();
        v_pdfs_floatParamsVaried.clear();
	int nit;
        if(par==0){
            for(int i=0; i<v_pdftype.size(); i++){
                if(_debug>=100)cout<<" vpdftype: "<<i<<"th --> "<<v_pdftype[i]<<endl;
                vrdm.push_back(-999);
                switch (v_pdftype[i]){
                    case typeLogNormal:
                    case typeShapeGaussianLinearMorph:
                    case typeShapeGaussianQuadraticMorph:
                        // for LHC-type frequentist method 
                        vrdm.back()=_rdm->Gaus(bUseBestEstimateToCalcQ!=2?0:(scaled?_fittedParsInData_sb[i]:_fittedParsInData_bonly[i]));
                        break;

                    case typeTruncatedGaussian:
                        //      another way is to throw normal gaus random number and regenerate if x<-1, it's more transparent
			tmp = _rdm->Gaus(bUseBestEstimateToCalcQ!=2?0:(scaled?_fittedParsInData_sb[i]:_fittedParsInData_bonly[i]), 1);
			nit = 0;
			while(tmp * v_TruncatedGaussian_maxUnc[i] < -1 ) {
				// FIXME  need to be care about range of nuisances, since we change the central values
				//tmp = _rdm->Gaus(bUseBestEstimateToCalcQ!=2?0:(scaled?_fittedParsInData_sb[i]:_fittedParsInData_bonly[i]), v_TruncatedGaussian_maxUnc[i]);
				tmp = _rdm->Gaus(bUseBestEstimateToCalcQ!=2?0:(scaled?_fittedParsInData_sb[i]:_fittedParsInData_bonly[i]), 1);
				nit++;
				if(nit > 1000 ) break;
			}
			vrdm.back()=tmp;
			break;

		    case typeGamma:
			//vrdm.back()=_rdm->Gamma(v_GammaN[i]);
			//frequentist way, need to toss integer number 
			if(_tossToyConvention)vrdm.back()=(int)_rdm->Gamma(bUseBestEstimateToCalcQ!=2?v_GammaN[i]:(scaled?_fittedParsInData_sb[i]:_fittedParsInData_bonly[i])); 
                        else vrdm.back()=_rdm->Gamma(bUseBestEstimateToCalcQ!=2?v_GammaN[i]:(scaled?_fittedParsInData_sb[i]:_fittedParsInData_bonly[i]));
                        if(_debug>=10)cout<<"DELETEME gamma N="<<(bUseBestEstimateToCalcQ!=2?v_GammaN[i]:(scaled?_fittedParsInData_sb[i]:_fittedParsInData_bonly[i]))<<endl;
                        if(_debug>=10)cout<<"DELETEME -->rdm = "<<vrdm.back()<<endl;
                        break;
                    case typeBifurcatedGaussian:
                        {
                            // FIXME need to think about this .. .
                            if(_debug>=100) cout<<" generating new value for parameter "<<v_uncname[i-1]<<endl;
                            tmprrv = _w_varied->var(TString::Format("%s_mean",v_uncname[i-1].c_str()));
                            //tmprrv->setDirtyInhibit(0);
                            tmprrv->setVal(bUseBestEstimateToCalcQ!=2?v_pdfs_floatParamsUnc[i][0]:(scaled?_fittedParsInData_sb[i]:_fittedParsInData_bonly[i]));

                            RooAbsData *tmpRDS = (RooAbsData*) _w_varied->pdf(TString::Format("%s_bfg",v_uncname[i-1].c_str()))
                                ->generate(RooArgSet(*_w_varied->var(TString::Format("%s_x", v_uncname[i-1].c_str()))), 1);
                            vrdm.back()= dynamic_cast<RooRealVar*> ( tmpRDS->get(0)->first() )->getVal();
                            delete tmpRDS;
                            if(_debug>=100) cout<<"  "<<vrdm.back()<<endl;
                            break;
                        }
                    case typeFlat:
                        //vrdm.back()=_rdm->Rndm(bUseBestEstimateToCalcQ!=2?0:(scaled?_fittedParsInData_sb[i]:_fittedParsInData_bonly[i]));
                        vrdm.back() = _rdm->Rndm(); // FIXME HGG  how to toss flat random number around fitted data
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
        }else{
            vrdm.push_back(0);
            for(int i=1; i<v_pdftype.size(); i++){
                vrdm.push_back(par[i]);
                if(_debug>=100)cout<<" index "<<i<<": "<<par[i]<<endl;
            }
        }		
        //if(_debug) cout<<"done for random gen"<<endl;

        VChannelVSample vv;
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
        if(includeCountingParts){
            if(bUseBestEstimateToCalcQ==1)vv= scaled?vv_exp_sigbkgs_scaled:vv_exp_sigbkgs;
            else if(bUseBestEstimateToCalcQ==0)vv= scaled?vv_randomized_sigbkgs_scaled:vv_randomized_sigbkgs;
            else vv= scaled?vv_fitted_sigbkgs_scaled:vv_fitted_sigbkgs;
	if(_debug >=1000 ){
		cout<<" DELETEME 88881 vv_fitted_sigbkgs.size()="<<vv_fitted_sigbkgs.size()<<endl;
		cout<<" DELETEME 88881 vv.size()="<<vv.size()<<endl;
        }
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
				if(_debug >=1000 ){
					cout<<"DELETEME ran = "<<ran<<" "<<endl;
					cout<<"DELETEME uncpars = "<<*uncpars<<" "<<endl;
					cout<<"DELETEME uncpars+1 = "<<*(uncpars+1)<<" "<<endl;
					cout<<"DELETEME -- "<<(ran>0? (*(uncpars+1)):(*uncpars))<<endl;
					cout<<"DELETEME -- "<<pow( (1+ (ran>0? (*(uncpars+1)):(*uncpars)) ) , ran )<<endl;
					cout<<" ch="<<ch<<" sam="<<isam<<endl; 
					cout<<"vv[ch][isam]="<<vv[ch][isam]<<endl;
				}
                                vv[ch][isam]*=pow( (1+ (ran>0? (*(uncpars+1)):(*uncpars)) ) , ran );
                                break;

                            case typeTruncatedGaussian :
                                //vv[ch][isam]*=( 1+vvvv_uncpar[ch][isam][iunc][(ran>0?1:0)] / v_TruncatedGaussian_maxUnc[indexcorrl] * ran );		
                                vv[ch][isam]*=( 1+vvvv_uncpar[ch][isam][iunc][(ran>0?1:0)] * ran );		
                                break;

                            case typeGamma :
                                if(vvvv_uncpar[ch][isam][iunc][0]>0){
                                    tmp = scaled?vv_exp_sigbkgs_scaled[ch][isam]:vv_exp_sigbkgs[ch][isam];
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
                            case typeFlat:
                                //cout<<"Rndm = "<<ran<<endl;
                                vv[ch][isam]*=( vvvv_uncpar[ch][isam][iunc][0] + (vvvv_uncpar[ch][isam][iunc][1]-vvvv_uncpar[ch][isam][iunc][0])*ran );
                                //cout<<"Rndm -> b = "<<vv[ch][isam]<<endl;
                                // 0: min,  1: max, 2: range
                                //cout<<"WARNING: typeFlat pdf on normalization not yet implemented, skip it "<<endl;
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
        }

        if(bHasParametricShape){
            //vv_pdfs_params_varied = bUseBestEstimateToCalcQ?vv_pdfs_params:vv_pdfs_params_randomized;
            if(bUseBestEstimateToCalcQ==1)vv_pdfs_norm_varied = scaled?vv_pdfs_norm_scaled:vv_pdfs_norm;
            else if(bUseBestEstimateToCalcQ==0)vv_pdfs_norm_varied = scaled?vv_pdfs_norm_randomized_scaled:vv_pdfs_norm_randomized;
            else vv_pdfs_norm_varied = scaled?vv_pdfs_norm_fitted_scaled:vv_pdfs_norm_fitted;
            for(int ch=0; ch<vv_pdfs.size(); ch++){
                nsigproc = v_pdfs_sigproc[ch];
                for(isam=0; isam<vv_pdfs[ch].size(); isam++){
                    for(iunc=0; iunc<vvv_pdfs_idcorrl[ch][isam].size(); iunc++){
                        indexcorrl = vvv_pdfs_idcorrl[ch][isam][iunc];
                        pdftype = vvv_pdfs_pdftype[ch][isam][iunc];
                        ran = vrdm[indexcorrl];
                        switch (pdftype){
                            case typeLogNormal:
                                if(_debug>=100) cout<<" in FluctuatedNumbers, bHasParametricShape: ch["<<ch<<"] isam["<<isam<<"] iunc["<<iunc<<"]"<<" pdftype="<<pdftype<<endl;
                                vv_pdfs_norm_varied[ch][isam] *= pow( (1+ (ran>0? vvv_pdfs_normvariation[ch][isam][iunc][1]: vvv_pdfs_normvariation[ch][isam][iunc][0])) , ran );
                                if(_debug>=100) cout<<" norm varied after this unc = "<<vv_pdfs_norm_varied[ch][isam]<<endl;
                                break;

                            case typeTruncatedGaussian :
                                if(_debug>=100) cout<<" in FluctuatedNumbers, bHasParametricShape: ch["<<ch<<"] isam["<<isam<<"] iunc["<<iunc<<"]"<<" pdftype="<<pdftype<<endl;
                                //vv_pdfs_norm_varied[ch][isam]*=( 1+vvv_pdfs_normvariation[ch][isam][iunc][(ran>0?1:0)] / v_TruncatedGaussian_maxUnc[indexcorrl] * ran );		
                                vv_pdfs_norm_varied[ch][isam]*=( 1+vvv_pdfs_normvariation[ch][isam][iunc][(ran>0?1:0)] * ran );		
                                if(_debug>=100) cout<<" norm varied after this unc = "<<vv_pdfs_norm_varied[ch][isam]<<endl;
                                break;

                            case typeGamma:
                                if(_debug>=100) cout<<" in FluctuatedNumbers, bHasParametricShape: ch["<<ch<<"] isam["<<isam<<"] iunc["<<iunc<<"]"<<" pdftype="<<pdftype<<endl;
                                if(vvv_pdfs_normvariation[ch][isam][iunc][0]>0){
                                    tmp = scaled?vv_pdfs_norm_scaled[ch][isam]:vv_pdfs_norm[ch][isam];
                                    if(isam<nsigproc){
                                        if(tmp!=0) {vv_pdfs_norm_varied[ch][isam] /=tmp; vv_pdfs_norm_varied[ch][isam]*=(ran * vvv_pdfs_normvariation[ch][isam][iunc][0] * _common_signal_strength );}
                                        else vv_pdfs_norm_varied[ch][isam] = ran * vvv_pdfs_normvariation[ch][isam][iunc][0] * _common_signal_strength ; // Gamma
                                    }else{
                                        if(tmp!=0) {vv_pdfs_norm_varied[ch][isam] /=tmp; vv_pdfs_norm_varied[ch][isam]*=(ran * vvv_pdfs_normvariation[ch][isam][iunc][0] );}
                                        else vv_pdfs_norm_varied[ch][isam] = ran * vvv_pdfs_normvariation[ch][isam][iunc][0]; // Gamma
                                    }
                                }else{ // if rho<0,   then this is multiplicative gamma function ....
                                    vv_pdfs_norm_varied[ch][isam] *= (ran/v_GammaN[indexcorrl]);
                                }
                                if(_debug>=100) cout<<" norm varied after this unc = "<<vv_pdfs_norm_varied[ch][isam]<<endl;
                                break;
                            case typeFlat:
                                cout<<"WARNING: typeFlat pdf on normalization not yet implemented, skip it "<<endl;
                                break;
                            default:
                                cout<<"ERROR: for shape norm uncertainies, we have only support LogNormal (lnN), TruncatedGaussian (trG) and Gamma (gmN) currently. exit "<<endl;
                                exit(1);
                                break;	
                        }
                    }
                    if(_debug>=100) cout<<" **norm varied after all unc = "<<vv_pdfs_norm_varied[ch][isam]<<endl;
                    tmprrv = _w_varied->var(vv_pdfs_normNAME[ch][isam]);
                    //tmprrv->setDirtyInhibit(1);
                    tmprrv->setVal(vv_pdfs_norm_varied[ch][isam]);
                }
            }

            int id;
            MapStrVV::iterator iter; 
            for(int i=0; i<v_pdfs_floatParamsName.size(); i++){

                int ip  = v_pdfs_floatParamsIndcorr[i];

                if(_debug>=100)cout<<v_pdfs_floatParamsName[i]<<endl;
                // additive implementation of multi-sources affecting params
                double param = vrdm[v_pdfs_floatParamsIndcorr[i]];

                if(v_pdfs_floatParamsType[i]==typeFlat){
                    // for flatParam: vector-->    0  norminal value, 1 max-min,  2 dumy, 3 min, 4 max 
                    //cout<<" flat "<<v_pdfs_floatParamsName[i]<<": "<<v_pdfs_floatParamsUnc[ip][0]<<" ["<<v_pdfs_floatParamsUnc[ip][3]<<","<<v_pdfs_floatParamsUnc[ip][4]<<"], "<< v_pdfs_floatParamsUnc[ip][1]<<", rdm= "<<param<<endl;
                    param =( v_pdfs_floatParamsUnc[ip][3] + param* v_pdfs_floatParamsUnc[ip][1] );
                    //cout<<" flat "<<v_pdfs_floatParamsUnc[ip][0]<<" --> "<<param<<endl;
                }


                if(_debug>=100)cout<<"DELETEME: param norminal value: "<<v_pdfs_floatParamsUnc[ip][0]<<endl;
                if(_debug>=100)cout<<"DELETEME: param only fitunc "<<param<<endl;

                iter = map_param_sources.find(v_pdfs_floatParamsName[i]);
                if (iter != map_param_sources.end() ) {
                    vector< vector<double> > vvtmp = iter->second;
                    for(int ii=0; ii<vvtmp.size(); ii++){
                        id = (int)vvtmp[ii][0];
                        if(v_pdftype[id]==typeBifurcatedGaussian) continue;
                        // FIXME HGG
                        param += vrdm[id]*(vrdm[id]<0?vvtmp[ii][1]:vvtmp[ii][2]);
                    }
                }else ; //std::cout << v_pdfs_floatParamsName[i]<<" is not in map_param_sources" << '\n';
                if(_debug>=100)cout<<"DELETEME: param after allunc: "<<param<<endl;

                // currently, we don't rethrow random numbers if accumulative param not in its range
                if(param<v_pdfs_floatParamsUnc[ip][3]) param=v_pdfs_floatParamsUnc[ip][3]; // ip  motivated by the code in AddUncertaintyOnShapeParam
                if(param>v_pdfs_floatParamsUnc[ip][4]) param=v_pdfs_floatParamsUnc[ip][4];

                if(_debug>=100)cout<<"DELETEME: param after range check: "<<param<<"  in ["<<v_pdfs_floatParamsUnc[ip][3]<<", "<<v_pdfs_floatParamsUnc[ip][4]<<"]"<<endl;
                if(_debug>=100) cout<<" setting new value for parameter "<<v_pdfs_floatParamsName[i]<<": "<<param<<endl;
                tmprrv = _w_varied->var(v_pdfs_floatParamsName[i].c_str());
                //tmprrv->setDirtyInhibit(1);
                tmprrv->setVal(param);
                v_pdfs_floatParamsVaried.push_back(param);

            }
            for(int ch=0; ch<vv_pdfs.size(); ch++){
                for(isam=0; isam<vv_pdfs[ch].size(); isam++){
                    double tmp = vv_pdfs_norm_varied[ch][isam];
                    if(_debug>=100) cout<<" **norm varied after all unc (only in datacard)= "<<tmp<<endl;
                    if(vv_pdfs_extranormNAME[ch][isam]!="") {
                        if(_debug>=100){ cout<<vv_pdfs_extranormNAME[ch][isam]<<endl;
                            if(static_cast<RooAbsReal*>(_w_varied->arg(vv_pdfs_extranormNAME[ch][isam]))) cout<<" exist "<<endl;
                            else cout<<" don't exist "<<endl;
                        }
                        if((RooAbsReal*)(_w_varied->arg(vv_pdfs_extranormNAME[ch][isam]))) {
                            tmp *= ((RooAbsReal*)(_w_varied->arg(vv_pdfs_extranormNAME[ch][isam])))->getVal();
                            if(_debug>=100) {
                                cout<<"      extranorm ="<< static_cast<RooAbsReal*>(_w_varied->arg(vv_pdfs_extranormNAME[ch][isam]))->getVal()<<endl;
                                cout<<" **norm varied after all unc (including extranorm)= "<<tmp<<endl;
                            }
                        }
                        vv_pdfs_norm_varied[ch][isam] = tmp;
                    }
                    tmprrv=_w_varied->var(vv_pdfs_normNAME[ch][isam]);
                    //tmprrv->setDirtyInhibit(1);
                    tmprrv->setVal(tmp);
                }
            }

            if(_debug>=100) {
                cout<<"FluctuatedNumbers, varied workspace: "<<endl;
                //_w_varied->Print("V"); // cause crash   ..... at some point
            }
        }

        if(par==0){
            //cout<<"DELETEME0 randomized parameters: "<<endl;
            for(int i=1; i<=max_uncorrelation; i++) {
                _randomizedPars[i] = vrdm[i]; // only for LHC-type test statistics evaluation
                //if(par==0)cout<<" par "<<i<<" "<<_randomizedPars[i]<<endl;
            }
        }
        return vv;
    }	
    VIChannel CountingModel::GetToyData_H0(double *pars){
        // background only hypothesis
        VChannelVSample vv = FluctuatedNumbers(pars);

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
                //v_pdfs_roodataset_toy.push_back( _w_varied->pdf(v_pdfs_b[ch])->generate(RooArgSet(*_w->var(v_pdfs_observables[ch])), Extended()));
		bool binned = false;
		if(v_pdfs_roodataset_real[ch]->InheritsFrom("RooDataHist")){binned = true;}
		if (binned)
			v_pdfs_roodataset_toy.push_back(_w_varied->pdf(v_pdfs_b[ch])->generateBinned(RooArgSet(*(v_pdfs_roodataset_real[ch]->get())), Extended()));
		else
			v_pdfs_roodataset_toy.push_back((RooAbsData*) _w_varied->pdf(v_pdfs_b[ch])->generate(RooArgSet(*(v_pdfs_roodataset_real[ch]->get())), Extended()));
                if(TString(v_pdfs_roodataset_toy[ch]->GetName())=="emptyData") {
                    if(_debug>=100)cout<<"changing toy name "<<endl;
                    v_pdfs_roodataset_toy[ch]->SetName( TString(v_pdfs_b[ch]+"Data").Data() );
                }

                if(_debug>=100){
                    v_pdfs_roodataset_toy[ch]->Print("V");
                    cout<<"toy name = "<<v_pdfs_roodataset_toy[ch]->GetName()<<endl;
                }

                if(_debug>=1000){
                    TString s = "H0_"; s+= (int)(10000000*_rdm->Rndm()); s+=".root";
                    TFile f(s, "RECREATE") ;
                    TH1F h("H0","H0", 100, _w->var(v_pdfs_observables[ch])->getMin(), _w->var(v_pdfs_observables[ch])->getMax() );
                    v_pdfs_roodataset_toy[ch]->fillHistogram(&h, RooArgList(*_w->var(v_pdfs_observables[ch])));
                    f.WriteTObject(&h);
                    f.Close();
                }

                if(_debug>=10)cout<<"H0, number of events generated for channel "<<ch<<": "<<v_pdfs_roodataset_toy[ch]->sumEntries()<<endl;
		/*
                for(int i=0; i<v_pdfs_roodataset_toy[ch]->sumEntries(); i++){
                    RooRealVar *r = dynamic_cast<RooRealVar*>(v_pdfs_roodataset_toy[ch]->get(i)->first());
                    vtmp.push_back( r->getVal() );
                }
                vv_pdfs_data_toy.push_back(vtmp);
		*/
            }
        }	

        return v;
    }
    VIChannel CountingModel::GetToyData_H1(double* pars){
        // alternative hypothesis
        VChannelVSample vv = FluctuatedNumbers(pars);

        if(_debug>=100)cout<<" in GetToyData_H1 ......"<<endl;

        VIChannel v; v.clear();
        double tmp;
        double br = GetSignalScaleFactor();
        for(int ch=0; ch<vv.size(); ch++){
            tmp=0;
            for(int isam=0; isam<vv[ch].size(); isam++){
                //start from 0, add signal
                tmp+=vv[ch][isam];	
            }

            if(_PhysicsModel == typeChargedHiggs){
                //****************** start ---  for charged higgs
                tmp=0;
                for(int isam=v_sigproc[ch]; isam<vv[ch].size(); isam++){
                    //start from 1,  don't add signal
                    tmp+=vv[ch][isam];	
                }

                // genelized for all channels
                // HH_ltau
                tmp +=  (br*br*vv[ch][0]); // HH -> r^2 * HH
                // WH_ltau
                tmp +=  (2*br*(1-br)*vv[ch][1]); //WH -> 2*r*(1-r)*WH
                // WW_ltau (tt->ltau, real tau) 
                tmp +=  ((1-br)*(1-br)*vv[ch][2]); //WW -> (1-r)*(1-r)*WW
                // WW_ltau (tt~->ll, also part of WW, one lepton fakes as tau) 
                tmp +=  ((1-br)*(1-br)*vv[ch][3]); //WW -> (1-r)*(1-r)*WW

                tmp -= vv[ch][2];
                tmp -= vv[ch][3];

                //********************** end --  for charged higgs
            }

            v.push_back(_rdm->Poisson(tmp));
        }

        if(_debug>=100)cout<<" in GetToyData_H1 done for counting channels"<<endl;
        if(bHasParametricShape){
            if(_debug>=100)cout<<" in GetToyData_H1 start to destroy old toy"<<endl;
            for(int ch=0; ch<v_pdfs_roodataset_toy.size(); ch++){
                //if(v_pdfs_roodataset_toy[ch])	delete v_pdfs_roodataset_toy[ch];
                delete v_pdfs_roodataset_toy[ch];
            }
            v_pdfs_roodataset_toy.clear();
            if(_debug>=100)cout<<" in GetToyData_H1 old toy destroyed"<<endl;
            for(int ch=0; ch<vv_pdfs.size(); ch++){
                vector<double> vtmp;vtmp.clear();

                //FIXME   we need release memory of used toys,  but delete them cause some segmentation fault 
                //if(!v_pdfs_roodataset_toy[ch]->IsZombie()) delete v_pdfs_roodataset_toy[ch];
                //v_pdfs_roodataset_toy[ch]->Delete();

                if(_debug>=100)cout<<" in GetToyData_H1, bHasParametricShape: ch = "<<ch<<endl;
                //v_pdfs_roodataset_toy.push_back(_w_varied->pdf(v_pdfs_sb[ch])->generate(RooArgSet(*_w->var(v_pdfs_observables[ch])), Extended()));
		bool binned = false;
		if(v_pdfs_roodataset_real[ch]->InheritsFrom("RooDataHist")){binned = true;}
		if (binned)
			v_pdfs_roodataset_toy.push_back(_w_varied->pdf(v_pdfs_sb[ch])->generateBinned(RooArgSet(*(v_pdfs_roodataset_real[ch]->get())), Extended()));
		else
			v_pdfs_roodataset_toy.push_back((RooAbsData*)_w_varied->pdf(v_pdfs_sb[ch])->generate(RooArgSet(*(v_pdfs_roodataset_real[ch]->get())), Extended()));
                // if generated events = 0, then the name of dataset will be emptyData, need to be changed
                if(TString(v_pdfs_roodataset_toy[ch]->GetName())=="emptyData") {
                    v_pdfs_roodataset_toy[ch]->SetName( TString(v_pdfs_sb[ch]+"Data").Data() );
                }

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
		/*
                for(int i=0; i<v_pdfs_roodataset_toy[ch]->sumEntries(); i++){
                    RooRealVar *r = dynamic_cast<RooRealVar*>(v_pdfs_roodataset_toy[ch]->get(i)->first());
                    vtmp.push_back( r->getVal() );
                }
                vv_pdfs_data_toy.push_back(vtmp);
		*/
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
        cout<<" [nofloat] nuisances: "<<endl;
        int nnn = 0;
        for(int i=0; i<v_uncname.size(); i++){
            if(v_uncFloatInFit[i]==false){
                cout<<i+1<<" "<<v_uncname[i]<<"  ";
                if(v_pdftype.size()>0) cout	<< v_pdftype[i+1] ;
                cout<<endl;
                nnn++;
            }
        }
        cout<< " totally "<<nnn<<" nuisances are not float in fit"<<endl;
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
    void CountingModel::SetSignalScaleFactor(double r, int bScaleBestEstimate){
        if(_debug>=10) cout<<"\n  *** SetSignalScaleFactor r= "<<r<<endl;
        if(!b_AllowNegativeSignalStrength && r<=0 ){
            cout<<"Error: signal strength r <=0"<<endl;
            cout<<"If you want to allow signal strength to be non-positive, please \n *** model->SetAllowNegativeSignalStrength(true)"<<endl;
            return;
        }
        _common_signal_strength=r;

        vv_exp_sigbkgs_scaled = vv_exp_sigbkgs;

        if(_PhysicsModel==typeChargedHiggs){
            // do nothing....   scaling signal will be done in somewhere else
        }else{
            for(int ch=0; ch<vv_exp_sigbkgs_scaled.size(); ch++){
                for(int isam=0; isam<v_sigproc[ch]; isam++){
                    vv_exp_sigbkgs_scaled[ch][isam]*=_common_signal_strength;
                }
            }
        }

        if(bHasParametricShape){
            RooRealVar*tmprrv;
            vv_pdfs_norm_scaled = vv_pdfs_norm;
            for(int ch=0; ch<vv_pdfs_norm_scaled.size(); ch++){
                for(int isam=0; isam<v_pdfs_sigproc[ch]; isam++){
                    vv_pdfs_norm_scaled[ch][isam]*=_common_signal_strength;

                    //FIXME HGG not sure what the following does ,  need to check and if nothing, then delete it
                    double tmp = vv_pdfs_norm_scaled[ch][isam];
                    if(vv_pdfs_extranormNAME[ch][isam]!="") {
                        tmp *= ((RooAbsReal*)(_w->arg(vv_pdfs_extranormNAME[ch][isam])))->getVal();
                    }
                    tmprrv=_w->var(vv_pdfs_normNAME[ch][isam]);
                    //tmprrv->setDirtyInhibit(1);
                    tmprrv->setVal(tmp);
                    tmprrv=_w_varied->var(vv_pdfs_normNAME[ch][isam]);
                    //tmprrv->setDirtyInhibit(1);
                    tmprrv->setVal(tmp);
                }
            }

        }

        //if(bHasParametricShape)SetTmpDataForUnbinned(v_pdfs_roodataset);// reset to data , need refit for mu=r  //FIXME  need move to somewhere else needed

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

	/*
	//if allow signal strength to be non-positive, then please make sure sig+bkgs >=0 in each channel 
	if(b_AllowNegativeSignalStrength){
	for(int ch=0; ch<vv_exp_sigbkgs_scaled.size(); ch++){
	double tot = 0;
	for(int p=0; p<vv_exp_sigbkgs_scaled[ch].size(); p++)	tot+=vv_exp_sigbkgs_scaled[ch][p];
	if(tot<0) {
	cout<<"Error: negative tot yield in channel "<<ch+1<<endl;
	cout<<"Please SetAllowNegativeSignalStrength(false)"<<endl;
	cout<<"Or we can think about setting a lower bound for the strength.."<<endl;
	return;
	}
	}
	}
	*/
    };
    VDChannel CountingModel::Get_AsimovData(int b){ // 0 for background only, 1 for s+b hypothesis
	    VDChannel vtmpdata = v_data; // could be pseudo-data for bands
	    if(b==0){
		    if(_debug)cout<<"		Get Asimov dataset: background only hypothesis"<<endl;
		    if(v_data_asimovb.size()>0) {
			    return v_data_asimovb;
		    }else{
			    for(int i=0; i<v_data.size(); i++){
				    vtmpdata[i]=0;
				    for(int b=v_sigproc[i]; b<vv_exp_sigbkgs.at(i).size(); b++){
					    vtmpdata[i]+= vv_exp_sigbkgs.at(i).at(b);
				    }
			    }
			    v_data_asimovb = vtmpdata;
		    }
	    }else if(b==1){
		    if(_debug)cout<<"		Get Asimov dataset: sig+bkg hypothesis"<<endl;
		    for(int i=0; i<v_data.size(); i++){
			    vtmpdata[i]=0;
			    for(int b=0; b<vv_exp_sigbkgs.at(i).size(); b++){
				    vtmpdata[i]+= vv_exp_sigbkgs.at(i).at(b);
			    }
		    }
	    }
	    return vtmpdata;
    }
    void CountingModel::ConstructAsimovData(int b, bool nominal){// for pdf shape channels
	    
	    if(b==0){
		    // FIXME, need to provide two options to construct asimov dataset,  nominal and fitted
		    if(!nominal) { // from fitted nuisances
			    double *pars = new double[cms_global->Get_max_uncorrelation()+1];
			    DoAfit(0, v_data_real, v_pdfs_roodataset_real, pars); 
			    //_fittedParsInData_bonly = pars;
			    Set_fittedParsInData_b(pars);

			    VChannelVSample vv = FluctuatedNumbers(_fittedParsInData_bonly);
			    v_data_asimovb = v_data;
			    for(unsigned int i=0; i<v_sigproc.size(); i++){
				    v_data_asimovb[i]=0;
				    for(unsigned int b=v_sigproc[i]; b<vv.at(i).size(); b++){
					    v_data_asimovb[i]+= vv.at(i).at(b);
				    }
			    }

			    int npars = cms_global->Get_max_uncorrelation();
			    for(int i=0; i<=npars; i++) 
				    printf("DELETEME  par %30s      %.6f \n", i>0?v_uncname[i-1].c_str():"signal_strength", _fittedParsInData_bonly[i]);
			    for(unsigned int i=0; i<v_data.size(); i++)
				    printf("DELETEME:  asimovb %d - %.5f\n", i, v_data_asimovb[i]);
		    }else{
				
			    v_data_asimovb = v_data;
			    for(unsigned int i=0; i<v_sigproc.size(); i++){
				    v_data_asimovb[i]=0;
				    for(unsigned int b=v_sigproc[i]; b<vv_exp_sigbkgs.at(i).size(); b++){
					    v_data_asimovb[i]+= vv_exp_sigbkgs.at(i).at(b);
				    }
			    }
		    }

		    v_pdfs_roodataset_asimovb.clear();
		    for(unsigned int i=0; i<v_pdfs_b.size(); i++){
			    //RooRealVar * observable = _w->var(v_pdfs_observables[i]);
			    RooArgSet * observable = (RooArgSet*) v_pdfs_roodataset_real[i]->get();
			    std::auto_ptr<TIterator> iter(observable->createIterator());
			    for(RooRealVar * obs= (RooRealVar*)iter->Next(); obs!=0; obs=(RooRealVar*)iter->Next()){
				    if(TString(obs->GetName()).BeginsWith("wgttmp_")) continue;
				    if(obs->getBins()>200) obs->setBins(200);
			    }
			    observable->Print();

			    RooDataHist *rdh_asimovb = _w_varied->pdf(v_pdfs_b[i]) -> generateBinned(*observable,ExpectedData());
			    //RooDataHist *rdh_asimovb = _w->pdf(v_pdfs_b[i]) -> generateBinned(*observable,ExpectedData());

			    TString swgt = "wgttmp_"; swgt+=v_pdfs_channelname[i];
			    //RooArgSet * rastmp = new RooArgSet(*observable, *(_w_varied->var(swgt)));
			    RooArgSet * rastmp = new RooArgSet(*observable, *(_w->var(swgt)));
			    RooDataSet *rds_asimovb = new RooDataSet(TString("rds_asimovb")+v_pdfs_channelname[i].c_str(), "rds_asimovb", *rastmp, swgt);
			    if(_debug)cout<<" rdh_asimovb["<<v_pdfs_channelname[i]<<"]->bins = "<<rdh_asimovb->numEntries()<<endl;
			    for(int j=0; j<rdh_asimovb->numEntries(); j++){

				    RooArgSet rastmp = *rdh_asimovb->get(j);
				    double wtmp = rdh_asimovb->weight();
				    rds_asimovb->add(rastmp, wtmp);
			    }
			    v_pdfs_roodataset_asimovb.push_back(rds_asimovb);


			    if(_debug>=10){
				    cout<<" *** DELETEME : "<<endl;
				    rds_asimovb->Print("v");
				    cout<<" sumEntries = "<<rds_asimovb->sumEntries()<<endl;
				    cout<<v_pdfs_channelname[i]<<"  asimovb numEntries = "<<rds_asimovb->numEntries()<<endl;
				    cout<<v_pdfs_channelname[i]<<"  asimovb sumEntries = "<<rds_asimovb->sumEntries()<<endl;

				    double xmin = _w->var(v_pdfs_observables[i])->getMin();
				    double xmax = _w->var(v_pdfs_observables[i])->getMax();
				    double nbins = _w->var(v_pdfs_observables[i])->getBinning().numBins();
				    double binwidth = (xmax-xmin)/(nbins);
				    TString s = "asimovb_"; s+= v_pdfs_channelname[i]; s+=".root";
				    TFile f(s, "RECREATE") ;
				    TH1F h("data","data", nbins, xmin, xmax);

				    for(int j=0; j<v_pdfs_roodataset_asimovb[i]->numEntries(); j++){
					    h.Fill(dynamic_cast<RooRealVar*>(v_pdfs_roodataset_asimovb[i]->get(j)->first())->getVal(), v_pdfs_roodataset_asimovb[i]->weight());
					    cout<<" c "<<i<<" bin "<<j<<": "<<dynamic_cast<RooRealVar*>(v_pdfs_roodataset_asimovb[i]->get(j)->first())->getVal()<<" weight "<<v_pdfs_roodataset_asimovb[i]->weight()<<endl;
				    }
				    v_pdfs_roodataset_asimovb[i]->Print("V");
				    v_pdfs_roodataset_asimovb[i]->Print("");
				    f.WriteTObject(&h);
			    }
		    }
	    }else if(b==1){
		    v_pdfs_roodataset_asimovsb.clear();
		    for(int i=0; i<v_pdfs_sb.size(); i++){
			    //RooRealVar * observable = _w->var(v_pdfs_observables[i]);
			    RooArgSet * observable = (RooArgSet*)v_pdfs_roodataset_real[i]->get();
			    std::auto_ptr<TIterator> iter(observable->createIterator());
			    for(RooRealVar * obs= (RooRealVar*)iter->Next(); obs!=0; obs=(RooRealVar*)iter->Next()){
				    if(TString(obs->GetName()).BeginsWith("wgttmp_")) continue;
				    if(obs->getBins()>200) obs->setBins(200);
				    cout<<"DELETEME ***** obs "<<obs->GetName()<<" bins = "<<obs->getBins()<<"\n";
			    }
			    RooDataHist *rdh_asimovsb = _w->pdf(v_pdfs_sb[i]) -> generateBinned(*observable,ExpectedData());
			    TString swgt = "wgttmp_"; swgt+=v_pdfs_channelname[i];
			    RooArgSet * rastmp = new RooArgSet(*observable, *(_w->var(swgt)));
			    RooDataSet *rds_asimovsb = new RooDataSet(TString("rds_asimovsb")+v_pdfs_channelname[i], "rds_asimovsb", *rastmp, swgt);
			    if(_debug)cout<<" rdh_asimovsb["<<v_pdfs_channelname[i]<<"]->bins = "<<rdh_asimovsb->numEntries()<<endl;
			    for(int j=0; j<rdh_asimovsb->numEntries(); j++){
				    RooArgSet rastmp = *rdh_asimovsb->get(j);
				    double wtmp = rdh_asimovsb->weight();
				    rds_asimovsb->add(rastmp, wtmp);
				    if(_debug && j==0)cout<<" DELETEME 3 "<<j<<" weight = "<<rdh_asimovsb->weight()<<endl;	
			    }
			    if(_debug>=10){
				    cout<<" *** DELETEME : "<<endl;
				    rds_asimovsb->Print("v");
				    cout<<" sumEntries = "<<rds_asimovsb->sumEntries()<<endl;
			    }
			    v_pdfs_roodataset_asimovsb.push_back(rds_asimovsb);
			    cout<<v_pdfs_channelname[i]<<"  asimovsb numEntries = "<<rds_asimovsb->numEntries()<<endl;
			    cout<<v_pdfs_channelname[i]<<"  asimovsb sumEntries = "<<rds_asimovsb->sumEntries()<<endl;
		    }
	    }
    }
    void CountingModel::UseAsimovData(int b){  // 0 for background only, 1 for s+b hypothesis
	    if(b==0){
		    cout<<"		Using Asimov dataset: background only hypothesis"<<endl;
		    v_data = Get_AsimovData(0);
		    v_data_real = Get_AsimovData(0);
		    v_pdfs_roodataset = v_pdfs_roodataset_asimovb;
		    v_pdfs_roodataset_tmp = v_pdfs_roodataset_asimovb;
	    }else if(b==1){
		    cout<<"		Using Asimov dataset: sig+bkg hypothesis"<<endl;
		    for(int i=0; i<v_data.size(); i++){
			    v_data[i]=0;
			    v_data_real[i]=0;
			    for(int b=0; b<vv_exp_sigbkgs.at(i).size(); b++){
				    v_data[i]+= vv_exp_sigbkgs.at(i).at(b);
				    v_data_real[i]+= vv_exp_sigbkgs.at(i).at(b);
			    }
		    }
		    v_pdfs_roodataset = v_pdfs_roodataset_asimovsb;
		    v_pdfs_roodataset_tmp = v_pdfs_roodataset_asimovsb;

	    }
    }
    CountingModel* CombineModels(CountingModel *cms1, CountingModel *cms2){
	    CountingModel *cms = new CountingModel();
	    if(cms1->GetDebug()>=10)cms->SetDebug(1);

	    VChannelVSampleVUncertaintyVParameter tmp_vvvv_uncpar = cms1->Get_vvvv_uncpar();
	    VChannelVSampleVUncertainty tmp_vvv_pdftype=cms1->Get_vvv_pdftype();	
	    VChannelVSampleVUncertainty tmp_vvv_idcorrl=cms1->Get_vvv_idcorrl();	
	    vector<string> tmp_v_uncname = cms1->Get_v_uncname();
	    vector<bool> tmp_v_uncFloatInFit = cms1->Get_v_uncFloatInFit();
	    if(cms2->GetDebug())cout<<"DELETEME 1"<<endl;
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
				    cms->TagUncertaintyFloatInFit(tmp_v_uncname[tmp_vvv_idcorrl.at(ch).at(isamp).at(iunc)-1], tmp_v_uncFloatInFit[tmp_vvv_idcorrl.at(ch).at(isamp).at(iunc)-1]);
			    }
		    }	
		    cms->AddObservedData(ch, cms1->Get_v_data()[ch]);
		    cms->SetProcessNames(ch, cms1->GetProcessNames(ch));
	    }	

	    tmp_vvvv_uncpar = cms2->Get_vvvv_uncpar();
	    tmp_vvv_pdftype=cms2->Get_vvv_pdftype();	
	    tmp_vvv_idcorrl=cms2->Get_vvv_idcorrl();	
	    tmp_v_uncname = cms2->Get_v_uncname();
	    tmp_v_uncFloatInFit = cms2->Get_v_uncFloatInFit();
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
				    cms->TagUncertaintyFloatInFit(tmp_v_uncname[tmp_vvv_idcorrl.at(ch).at(isamp).at(iunc)-1], tmp_v_uncFloatInFit[tmp_vvv_idcorrl.at(ch).at(isamp).at(iunc)-1]);
			    }
		    }	
		    cms->AddObservedData(newch, cms2->Get_v_data()[ch]);
		    cms->SetProcessNames(newch, cms2->GetProcessNames(ch));
	    }	

	    if(cms2->GetDebug()) cout<<"Added Counting parts"<<endl;

	    vector<CountingModel*> ms; ms.push_back(cms1); ms.push_back(cms2);
	    for(int m=0; m<ms.size(); m++){
		    if(ms[m]->hasParametricShape()){
			    if(cms2->GetDebug())cout<<"DELETEME 2"<<endl;
			    RooWorkspace * w = ms[m]->GetWorkSpace();
			    vector<int> vsigproc = ms[m]->Get_v_pdfs_sigproc();
			    vector< vector<string> > vvpdfs = ms[m]->Get_vv_pdfs();
			    vector<TString> vobs = ms[m]->Get_v_pdfs_observables();
			    vector< vector<double> > vnorm = ms[m]->Get_vv_pdfs_norm_nonscaled(); 
			    vector<string> vchnames = ms[m]->Get_v_pdfs_channelname();
			    vector< RooAbsData* > vrds = ms[m]->Get_v_pdfs_roodataset();
			    tmp_vvvv_uncpar = ms[m]->Get_vvv_pdfs_normvariation();
			    tmp_vvv_pdftype=ms[m]->Get_vvv_pdfs_pdftype();	
			    tmp_vvv_idcorrl=ms[m]->Get_vvv_pdfs_idcorrl();	
			    tmp_v_uncname = ms[m]->Get_v_uncname();
			    tmp_v_uncFloatInFit = ms[m]->Get_v_uncFloatInFit();
			    vector< vector<double> > vparamunc = ms[m]->Get_v_pdfs_floatParamsUnc(); // from 0 to max_uncorl
			    vector<int> vparamIndcorr = ms[m]->Get_v_pdfs_floatParamsIndcorr();      // only for params
			    vector<int> tmp_v_pdftype = ms[m]->Get_v_pdftype();
			    vector< vector<TString> > vvextranormname = ms[m]->Get_vv_pdfs_extranormNAME();
			    for(int ch=0; ch<vsigproc.size(); ch++){
				    int newch = cms->Get_vv_pdfs().size();
				    if(ms[m]->GetDebug()) cout<<"Adding ch = "<<newch<<"th channel from "<<ms[m]->GetModelName()<<endl;
				    vector<RooAbsPdf*> vs, vb;
				    vector<double> vsnorm, vbnorm;
				    vector<RooAbsArg*> vsExtraNorm, vbExtraNorm; 
				    for(int p =0 ; p<vvpdfs[ch].size(); p++){
					    if(p<vsigproc[ch]){
						    vs.push_back(w->pdf(vvpdfs[ch][p].c_str()));
						    vsnorm.push_back(vnorm[ch][p]);
						    vsExtraNorm.push_back((RooAbsArg*)w->arg(vvextranormname[ch][p]));
					    }
					    else{
						    vb.push_back(w->pdf(vvpdfs[ch][p].c_str()));
						    vbnorm.push_back(vnorm[ch][p]);
						    vbExtraNorm.push_back((RooAbsArg*)w->arg(vvextranormname[ch][p]));
					    }
				    }
				    RooRealVar*x=w->var(vobs[ch]);
				    cms->AddChannel(vchnames[ch], x, vs, vsnorm,vsExtraNorm, vb, vbnorm, vbExtraNorm);
				    if(cms2->GetDebug()) cout<<"  AddChannel "<<endl;
				    cms->AddObservedDataSet(vchnames[ch], vrds[ch]);
				    if(cms2->GetDebug()) cout<<"  AddObservedDataSet"<<endl;

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
						    }else if(tmp_vvv_pdftype.at(ch).at(isamp).at(iunc)==typeControlSampleInferredLogNormal
								    || tmp_vvv_pdftype.at(ch).at(isamp).at(iunc)==typeGamma){
							    cms->AddUncertaintyOnShapeNorm(vchnames[ch], isamp, 
									    tmp_vvvv_uncpar.at(ch).at(isamp).at(iunc).at(0), 
									    tmp_vvvv_uncpar.at(ch).at(isamp).at(iunc).at(1), 
									    tmp_vvvv_uncpar.at(ch).at(isamp).at(iunc).at(2), 
									    tmp_vvv_pdftype.at(ch).at(isamp).at(iunc),
									    tmp_v_uncname[tmp_vvv_idcorrl.at(ch).at(isamp).at(iunc)-1]
									    );

						    }else {
							    cout<<"The pdftype Not implemented yet in shapeNorm: "<<tmp_vvv_pdftype.at(ch).at(isamp).at(iunc)<<endl;
						    }
						    cms->TagUncertaintyFloatInFit(tmp_v_uncname[tmp_vvv_idcorrl.at(ch).at(isamp).at(iunc)-1], tmp_v_uncFloatInFit[tmp_vvv_idcorrl.at(ch).at(isamp).at(iunc)-1]);
					    }
				    }
				    if(cms2->GetDebug()) cout<<"  AddUncertaintyOnShapeNorm"<<endl;
				    if(cms2->GetDebug()) cout<<"Added ch = "<<newch<<"th channel from "<<ms[m]->GetModelName()<<endl;
			    }

			    // add uncertainties on shape parameters
			    if(cms2->GetDebug()) cout<<" uncname.size="<<tmp_v_uncname.size()<<" vparamunc.size="<<vparamunc.size()<<endl;
			    if(cms2->GetDebug()) cout<<" "<<vparamIndcorr.size()<<" floating parameters "<<endl;
			    for(int ii=0; ii<vparamIndcorr.size(); ii++){
				    if(cms2->GetDebug()>=10) cout<<" ii =  "<<ii<<endl;
				    int id= vparamIndcorr[ii];
				    if(cms2->GetDebug()>=10) cout<<" id =  "<<id<<endl;

				    if(tmp_v_pdftype[id]==typeBifurcatedGaussian)cms->AddUncertaintyOnShapeParam(tmp_v_uncname[id-1], vparamunc[id][0], vparamunc[id][1], vparamunc[id][2], vparamunc[id][3], vparamunc[id][4] );
				    else if (tmp_v_pdftype[id]==typeFlat) cms->AddUncertaintyOnShapeParam(tmp_v_uncname[id-1]);

				    if(cms2->GetDebug()) cout<<"  AddUncertaintyOnShapeParam "<<tmp_v_uncname[id-1]<<" "
					    <<vparamunc[id][0] <<" "
						    <<vparamunc[id][1] <<" "
						    <<vparamunc[id][2] <<" "
						    <<vparamunc[id][3] <<" "
						    <<vparamunc[id][4] <<" "
						    <<endl;
				    cms->TagUncertaintyFloatInFit(tmp_v_uncname[id-1], tmp_v_uncFloatInFit[id-1]);
			    }

			    if(cms2->GetDebug()>=10) cout<<" before adding map_param_sources"<<endl;

			    // add map_param_sources
			    MapStrVV map = ms[m]->Get_map_param_sources();

			    if(cms2->GetDebug()>=10) cout<<" map_param_sources.size = "<<map.size()<<endl;

			    MapStrVV::iterator iter = map.begin();
			    for(; iter!=map.end(); iter++){
				    if(cms2->GetDebug()>=10) cout<<" pname = "<<iter->first<<endl;
				    vector< vector<double> > vv = iter->second;
				    if(ms[m]->GetDebug()>=10) cout<<" vv.size = "<<vv.size()<<endl;
				    for(int ii=0; ii<vv.size(); ii++){
					    int id = (int)vv[ii][0];
					    cms->AddUncertaintyAffectingShapeParam(tmp_v_uncname[id-1], iter->first, vv[ii][1], vv[ii][2]);
					    cms->TagUncertaintyFloatInFit(tmp_v_uncname[id-1], tmp_v_uncFloatInFit[id-1]);
				    }
			    }
			    if(cms2->GetDebug()>=10) cout<<" after adding map_param_sources"<<endl;
		    }
	    }



	    return cms;

    }
    double CountingModel::GetExpectedNumber(int index_channel, int index_sample){
	    if(index_channel >= vv_exp_sigbkgs_scaled.size()) return -1;
	    if(index_sample>= vv_exp_sigbkgs_scaled.at(index_channel).size()) return -1;
	    return vv_exp_sigbkgs_scaled.at(index_channel).at(index_sample);
    }
    int CountingModel::RemoveChannelsWithExpectedSignal0orBkg0(int king){

	//	return 0;

	    // should be advoked at the end of model construction, otherwise you will get into trouble ....
	    // either before or after "ConfigUncertaintyPdfs"
	    std::vector< vector<double> >::iterator iter=vv_exp_sigbkgs.begin(); 
	    int skippedchannels_s = 0, skippedchannels_b=0, skippedchannels_sb=0;
	    int skippedchannels =0;
	    for(; iter!=vv_exp_sigbkgs.end();){
		    int position=iter-vv_exp_sigbkgs.begin();
		    double bkg = 0, sig=0;
		    for(int p=0; p<(*iter).size(); p++){
			    if( p>=v_sigproc[position] )bkg += (*iter)[p];
			    else sig += (*iter)[p];
		    }
		    if( 
				    ( (king==1 || king==2) && sig<=0 ) ||
				    ( (king==0 || king==2) && bkg<=0        ) 
		      ) {
			    iter=vv_exp_sigbkgs.erase( iter );
			    v_data.erase( v_data.begin()+position );
			    v_data_real.erase( v_data_real.begin()+position );
			    v_sigproc.erase( v_sigproc.begin()+position );
			    v_channelname.erase(v_channelname.begin()+position);
			    vv_exp_sigbkgs_scaled.erase( vv_exp_sigbkgs_scaled.begin()+position );
			    vvvv_uncpar.erase( vvvv_uncpar.begin()+position );
			    vvv_pdftype.erase( vvv_pdftype.begin()+position );
			    vvv_idcorrl.erase( vvv_idcorrl.begin()+position );
			    vv_procname.erase( vv_procname.begin()+position );
			    vvv_shapeuncindex.erase(vvv_shapeuncindex.begin()+position );

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
	    return skippedchannels;
    }


    // for parametric shapes
    void CountingModel::AddChannel(string channel_name, RooRealVar* observable, vector<RooAbsPdf*> sigPdfs, vector<double> sigNorms, vector<RooAbsArg*> vsExtraNorm,
		    vector<RooAbsPdf*> bkgPdfs, vector<double> bkgNorms, vector<RooAbsArg*> vbExtraNorm){
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
	    _w->import(*observable, RecycleConflictNodes());
	    _w_varied->import(*observable, RecycleConflictNodes());
	    v_pdfs_observables.push_back(observable->GetName());

	    if(_w->var(v_pdfs_observables.back())->getBins()>200){
		    _w->var(v_pdfs_observables.back())->setBins(200);
		    _w_varied->var(v_pdfs_observables.back())->setBins(200);
		    observable->setBins(200);
	    }


	    vector<string> vsigbkgs; vsigbkgs.clear();
	    vector< vector< vector<double> > > vvvuncpar; vvvuncpar.clear();
	    vector< vector<int> > vvpdftype; vvpdftype.clear();
	    vector< vector<int> > vvidcorrl; vvidcorrl.clear();
	    vector<double> vnorms; vnorms.clear();

	    vector< vector<double> > vvunc; vvunc.clear();
	    vector<int> vpdftype; vpdftype.clear();
	    vector<int> vidcorrl; vidcorrl.clear();

	    vector<TString> vrrvnorm; vrrvnorm.clear();
	    vector<TString> vextranorm; vextranorm.clear();

	    vector<TString> vcoeffs; vcoeffs.clear();

	    vector<RooRealVar*> vrrvparams; vrrvparams.clear();
	    vector< vector<RooRealVar*> > vvrrvparams; vvrrvparams.clear();
	    for(int i=0; i<sigPdfs.size(); i++){
		    vsigbkgs.push_back(sigPdfs[i]->GetName());
		    vvvuncpar.push_back(vvunc);
		    vvpdftype.push_back(vpdftype);
		    vvidcorrl.push_back(vidcorrl);

		    vnorms.push_back(sigNorms[i]);
		    TString sn = channel_name; sn+=vproc[i]; sn+="_norm";sn+=i;
		    RooRealVar *rrv = new RooRealVar(sn, "", sigNorms[i]);//FIXME HGG
		    vrrvnorm.push_back(sn);


		    if(vsExtraNorm[i]) { 
			    _w->import(*vsExtraNorm[i], RecycleConflictNodes());
			    _w_varied->import(*vsExtraNorm[i], RecycleConflictNodes());
			    vextranorm.push_back(vsExtraNorm[i]->GetName());
			    //cout<<vextranorm.back()<<endl;
			    if((RooAbsReal*)(_w_varied->arg(vextranorm.back()))) {
				    double tmp = ((RooAbsReal*)(_w_varied->arg(vextranorm.back())))->getVal();
				    if(_debug>=100) cout<<" **extranorm= "<<tmp<<endl;
			    }else{
				    cout<<" don't exist "<<endl;
			    }
			    rrv->setVal(sigNorms[i]*(((RooAbsReal*)(vsExtraNorm[i]))->getVal()));
		    }else{
			    vextranorm.push_back("");
		    }

		    _w->import(*rrv);
		    _w->import(*sigPdfs[i], RenameConflictNodes(channel_name.c_str()));

		    _w_varied->import(*rrv);
		    _w_varied->import(*sigPdfs[i], RenameConflictNodes(channel_name.c_str()));
		    //RooArgSet *rds	= sigPdfs[i]->getParameters(*observable);
		    // need to store the list of parameters and for future modification, fluctuation 


		    TString coef;
		    coef = rrv->GetName();
		    vcoeffs.push_back(coef);

		    if(_debug>=10){
			    if(bRedefineObservableRange){
				    _w->var(v_pdfs_observables.back())->setMin(ObservableRangeMin);
				    _w->var(v_pdfs_observables.back())->setMax(ObservableRangeMax);
				    _w->var(v_pdfs_observables.back())->setBins(ObservableBins);
			    }
			    TString s = "pdf_"; s+= channel_name; s+="_signals.root";
			    TFile f(s, i==0?"RECREATE":"UPDATE") ;
			    double xmin = _w->var(v_pdfs_observables.back())->getMin();
			    double xmax = _w->var(v_pdfs_observables.back())->getMax();
			    TString sp = channel_name; sp+=sigPdfs[i]->GetName();
			    int nbins = 100; if ( bRedefineObservableRange ) nbins = ObservableBins;
			    TH1F hs(sp,sp, nbins, xmin, xmax);
			    RooArgSet vars(*(_w->var(v_pdfs_observables.back() ) ) );
			    sigPdfs[i]->fillHistogram(&hs, vars);
			    hs.SetNormFactor(rrv->getVal());
			    hs.Write();
			    f.Close();
			    hs.Write();
			    f.Close();
		    }

	    }

	    for(int i=0; i<bkgPdfs.size(); i++){
		    vsigbkgs.push_back(bkgPdfs[i]->GetName());
		    vvvuncpar.push_back(vvunc);
		    vvpdftype.push_back(vpdftype);
		    vvidcorrl.push_back(vidcorrl);

		    tmp_totbkg+=bkgNorms[i];
		    vnorms.push_back(bkgNorms[i]);

		    TString sn = channel_name; sn+=vproc[i+sigNorms.size()]; sn+="_norm";sn+=(i+sigNorms.size());
		    RooRealVar *rrv = new RooRealVar(sn, "", bkgNorms[i]); // FIXME HGG
		    vrrvnorm.push_back(sn);

		    if(vbExtraNorm[i]) { 
			    _w->import(*vbExtraNorm[i]);
			    _w_varied->import(*vbExtraNorm[i]);
			    vextranorm.push_back(vbExtraNorm[i]->GetName());
			    rrv->setVal(bkgNorms[i]*(((RooAbsReal*)vbExtraNorm[i])->getVal()));
		    }else{
			    vextranorm.push_back("");
		    }

		    _w->import(*rrv);
		    _w->import(*bkgPdfs[i], RenameConflictNodes(channel_name.c_str()));

		    _w_varied->import(*rrv);
		    _w_varied->import(*bkgPdfs[i], RenameConflictNodes(channel_name.c_str()));

		    TString coef;
		    coef = rrv->GetName();
		    vcoeffs.push_back(coef);

		    if(_debug>=10){
			    if(bRedefineObservableRange){
				    _w->var(v_pdfs_observables.back())->setMin(ObservableRangeMin);
				    _w->var(v_pdfs_observables.back())->setMax(ObservableRangeMax);
				    _w->var(v_pdfs_observables.back())->setBins(ObservableBins);
			    }
			    TString s = "pdf_"; s+= channel_name; s+="_bkgs.root";
			    TFile f(s, i==0?"RECREATE":"UPDATE") ;
			    double xmin = _w->var(v_pdfs_observables.back())->getMin();
			    double xmax = _w->var(v_pdfs_observables.back())->getMax();
			    TString sp = channel_name; sp+=bkgPdfs[i]->GetName();
			    int nbins = 100; if ( bRedefineObservableRange ) nbins = ObservableBins;
			    TH1F hs(sp,sp, nbins, xmin, xmax);
			    RooArgSet vars(*(_w->var(v_pdfs_observables.back() ) ) );
			    bkgPdfs[i]->fillHistogram(&hs, vars);
			    hs.SetNormFactor(rrv->getVal());
			    hs.Write();
			    f.Close();
		    }

	    }

	    vv_pdfs.push_back(vsigbkgs);
	    vvv_pdfs_normvariation.push_back(vvvuncpar);
	    vvv_pdfs_pdftype.push_back(vvpdftype);
	    vvv_pdfs_idcorrl.push_back(vvidcorrl);
	    vv_pdfs_normNAME.push_back(vrrvnorm);
	    vv_pdfs_extranormNAME.push_back(vextranorm);
	    vv_pdfs_norm.push_back(vnorms);
	    vv_pdfs_norm_scaled.push_back(vnorms);
	    vv_pdfs_norm_varied.push_back(vnorms);
	    if(_debug>=100)cout<<"channel "<<channel_name<<": tot.size="<<signal_processes+bkg_processes<<" bkg.size="<<bkg_processes<<endl;


	    // construct  model_sb,  model_s,  model_b
	    TString s = "SUM::"; s+=channel_name; s+="_sb(";
	    for(int i=0; i<vsigbkgs.size(); i++){
		    if(i!=0) s+=",";
		    s+=vcoeffs[i];
		    s+="*"; s+=vsigbkgs[i]; 			
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
		    s+=vcoeffs[i];
		    s+="*"; s+=vsigbkgs[i]; 			
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
		    s+=vcoeffs[i];
		    s+="*"; s+=vsigbkgs[i]; 			
	    }
	    s+=")";
	    _w->factory(s);
	    _w_varied->factory(s);
	    s = channel_name; s+="_b";
	    v_pdfs_b.push_back(s);	
	    //if(_debug) _w->pdf(s)->Print("V");

	    TString swgt = "wgttmp_"; swgt+=channel_name;
	    if(_w->var(swgt) == NULL)_w->factory(swgt+"[0,10000000]");
	    //RooAbsData * rds = _w->pdf(v_pdfs_b.back()) -> generate(*observable, Extended());
	    RooAbsData * rds = (RooAbsData*)_w->pdf(v_pdfs_b.back()) -> generate(*observable, 1) ;
	    v_pdfs_roodataset.push_back(rds);
	    v_pdfs_roodataset_tmp.push_back(rds);
	    v_pdfs_roodataset_real.push_back(rds);
	    /*
	       RooDataHist *rdh_asimovb = _w->pdf(v_pdfs_b.back()) -> generateBinned(*observable,ExpectedData());
	       RooDataHist *rdh_asimovsb = _w->pdf(v_pdfs_sb.back()) -> generateBinned(*observable,ExpectedData());

	       RooArgSet * rastmp = new RooArgSet(*observable, *(_w->var("wgttmp")));
	       RooAbsData *rds_asimovb = new RooAbsData(TString("rds_asimovb")+channel_name.c_str(), "rds_asimovb", *rastmp, "wgttmp");
	       if(_debug)cout<<" rdh_asimovb["<<channel_name<<"]->bins = "<<rdh_asimovb->numEntries()<<endl;
	       for(int i=0; i<rdh_asimovb->numEntries(); i++){
	       if(_debug && i==0)cout<<" DELETEME 3 "<<i<<" weight = "<<rdh_asimovb->weight()<<endl;	
	       rds_asimovb->add(*rdh_asimovb->get(i), rdh_asimovb->weight());
	       }
	       if(_debug>=10){
	       cout<<" *** DELETEME : "<<endl;
	       rds_asimovb->Print("v");
	       cout<<" sumEntries = "<<rds_asimovb->sumEntries()<<endl;
	       }
	       v_pdfs_roodataset_asimovb.push_back(rds_asimovb);

	       RooAbsData *rds_asimovsb = new RooAbsData(TString("rds_asimovsb")+channel_name.c_str(), "rds_asimovsb", *rastmp, "wgttmp");
	       if(_debug)cout<<" rdh_asimovsb["<<channel_name<<"]->bins = "<<rdh_asimovsb->numEntries()<<endl;
	       for(int i=0; i<rdh_asimovsb->numEntries(); i++){
	       if(_debug && i==0)cout<<" DELETEME 3 "<<i<<" weight = "<<rdh_asimovsb->weight()<<endl;	
	       rds_asimovsb->add(*rdh_asimovsb->get(i), rdh_asimovsb->weight());
	       }
	       if(_debug>=10){
	       cout<<" *** DELETEME : "<<endl;
	       rds_asimovsb->Print("v");
	       cout<<" sumEntries = "<<rds_asimovsb->sumEntries()<<endl;
	       }
	       v_pdfs_roodataset_asimovsb.push_back(rds_asimovsb);
	       */
	    /*
	       if(_debug>=10){
	       double dtmp = 0;
	       for(int i=0; i<rds_asimovb->numEntries(); i++){
	       rds_asimovb->get(i); 
	       dtmp+=rds_asimovb->weight();
	       }
	       cout<<" from numEntries and weight,  sumEntries  = "<<dtmp<<endl;
	       }
	       */
	    //v_pdfs_roodataset_toy.push_back(rds);
	    //cout<<"H0, number of events generated for channel "<<ch<<": "<<v_pdfs_roodataset_toy[ch]->sumEntries()<<endl;
	    if(_debug) 	rds->Print("V");
	    /*
	       vector<double> vdata;
	       for(int i=0; i<rds->sumEntries(); i++){
	       RooRealVar *r = dynamic_cast<RooRealVar*>(rds->get(i)->first());
	       vdata.push_back( r->getVal() );
	       }
	       vv_pdfs_data.push_back(vdata);
	       */

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
	    double weight = 1;

	    double btot = 0, stot=0;
	    for(int i=0; i<vv_pdfs_norm_scaled[ch].size(); i++){  //FIXME HGG   need to multiply the extra norm in the _w 
		    if(i>=v_pdfs_sigproc[ch]) btot+=vv_pdfs_norm_scaled[ch][i];
		    else stot+=vv_pdfs_norm_scaled[ch][i]; 
	    }
	    //RooArgSet vars(*(_w->var(v_pdfs_observables[ch])));
	    RooArgSet vars(*(v_pdfs_roodataset_real[ch])->get());
	    if(dataOrToy == 0){
		    int ntot = int(v_pdfs_roodataset[ch]->sumEntries());
		    for(int i=0; i<v_pdfs_roodataset[ch]->numEntries(); i++){
			   // _w->var(v_pdfs_observables[ch])->setVal(( dynamic_cast<RooRealVar*>(v_pdfs_roodataset[ch]->get(i)->first()))->getVal());
			    
			    std::auto_ptr<TIterator> iter(v_pdfs_roodataset[ch]->get(i)->createIterator());
			    for(RooRealVar * obs= (RooRealVar*)iter->Next(); obs!=0; obs=(RooRealVar*)iter->Next()){
				    if(TString(obs->GetName()).BeginsWith("wgttmp_")) continue;
				    _w->var(obs->GetName())->setVal(obs->getVal());
			    }

			    if(_debug>=100){
				    if(i==0 or i==ntot-1 or i==ntot/2){
					    cout<<"* event "<<i<<":  m= "<<( dynamic_cast<RooRealVar*>(v_pdfs_roodataset[ch]->get(i)->first()))->getVal()<<endl;
					    cout<<" pdfs= "<<_w->pdf(v_pdfs_s[ch])->getVal(&vars)<<endl;
					    cout<<" pdfb= "<<_w->pdf(v_pdfs_b[ch])->getVal(&vars)<<endl;
					    cout<<"stot = "<<stot<<" btot="<<btot<<endl;
				    }
			    }
			    weight = v_pdfs_roodataset[ch]->weight();
			    ret+= log(1+ stot/btot*_w->pdf(v_pdfs_s[ch])->getVal(&vars)/_w->pdf(v_pdfs_b[ch])->getVal(&vars)) * weight;
		    }
	    }else if(dataOrToy==1){
		    int ntot = int(v_pdfs_roodataset_toy[ch]->sumEntries());
		    for(int i=0; i<v_pdfs_roodataset_toy[ch]->numEntries(); i++){
			    ///_w->var(v_pdfs_observables[ch])->setVal(( dynamic_cast<RooRealVar*>(v_pdfs_roodataset_toy[ch]->get(i)->first()))->getVal());
			    
			    std::auto_ptr<TIterator> iter(v_pdfs_roodataset_toy[ch]->get(i)->createIterator());
			    for(RooRealVar * obs= (RooRealVar*)iter->Next(); obs!=0; obs=(RooRealVar*)iter->Next()){
				    if(TString(obs->GetName()).BeginsWith("wgttmp_")) continue;
				    _w->var(obs->GetName())->setVal(obs->getVal());
			    }
			    if(_debug>=100){
				    if(i==0 or i==ntot-1 or i==ntot/2){
					    cout<<"* event "<<i<<":  m= "<<( dynamic_cast<RooRealVar*>(v_pdfs_roodataset_toy[ch]->get(i)->first()))->getVal()<<endl;
					    cout<<" pdfs= "<<_w->pdf(v_pdfs_s[ch])->getVal(&vars)<<endl;
					    cout<<" pdfb= "<<_w->pdf(v_pdfs_b[ch])->getVal(&vars)<<endl;
					    cout<<"stot = "<<stot<<" btot="<<btot<<endl;
				    }
			    }
			    weight = v_pdfs_roodataset[ch]->weight();
			    ret+= log(1+ stot/btot*_w->pdf(v_pdfs_s[ch])->getVal(&vars)/_w->pdf(v_pdfs_b[ch])->getVal(&vars)) * weight;
		    }
	    }

	    ret-=stot;
	    if(_debug>=10)cout<<"EvaluateLnQ of "<<(dataOrToy==0?"data":"toy")<<" in channel ["<<v_pdfs_channelname[ch]<<"]: lnQ= "<<ret<<endl;
	    return ret;
    }

    void CountingModel::AddObservedDataSet(int index_channel, RooAbsData* rds){
	    //if(v_pdfs_roodataset[index_channel]) delete v_pdfs_roodataset[index_channel];
	    int ch = index_channel;

	    RooAbsData *tmp = v_pdfs_roodataset[ch];
	    v_pdfs_roodataset[ch]=rds;
	    v_pdfs_roodataset_tmp[ch]=rds;
	    v_pdfs_roodataset_real[ch]=rds;
	    cout<<v_pdfs_channelname[ch]<<" data numEntries = "<<v_pdfs_roodataset_real[ch]->numEntries()<<endl;;
	    cout<<v_pdfs_channelname[ch]<<" data sumEntries = "<<v_pdfs_roodataset_real[ch]->sumEntries()<<endl;
	    delete tmp;

	    if(_debug>=10){
		    TString s = "data_"; s+= v_pdfs_channelname[ch]; s+=".root";
		    TFile f(s, "RECREATE") ;
		    double xmin = _w->var(v_pdfs_observables[ch])->getMin();
		    double xmax = _w->var(v_pdfs_observables[ch])->getMax();
		    double nbins = _w->var(v_pdfs_observables[ch])->getBinning().numBins();
		    cout<<"  "<<v_pdfs_observables[ch]<<" binning "<<nbins<<endl;

		    TH1F hs("s","s", nbins, xmin, xmax);
		    TH1F hb("b","b", nbins,  xmin, xmax);
		    TH1F hsb("sb","sb", nbins,  xmin, xmax);
		    TH1F h("data","data", nbins, xmin, xmax);
		    double binwidth = (xmax-xmin)/nbins;
		    RooArgSet vars(*(_w->var(v_pdfs_observables.back() ) ) );
		    for(int x = 1; x<=nbins; x++){
			    _w->var(v_pdfs_observables.back())->setVal((xmax-xmin)/nbins*(x-1)+xmin); 
			    hs.SetBinContent(x, _w->pdf(v_pdfs_s.back())->getVal(&vars));
			    hb.SetBinContent(x, _w->pdf(v_pdfs_b.back())->getVal(&vars));
			    hsb.SetBinContent(x, _w->pdf(v_pdfs_sb.back())->getVal(&vars));
		    }

		    hs.Write();
		    hb.Write();
		    hsb.Write();

		    for(int i=0; i<v_pdfs_roodataset[ch]->numEntries(); i++){
			    //_w->var(v_pdfs_observables[ch])->setVal(( dynamic_cast<RooRealVar*>(v_pdfs_roodataset[ch]->get(i)->first()))->getVal());
			    //h.Fill(_w->var(v_pdfs_observables[ch])->getVal());
			    h.Fill(dynamic_cast<RooRealVar*>(v_pdfs_roodataset[ch]->get(i)->first())->getVal(), v_pdfs_roodataset[ch]->weight());
			    if(_debug>=10)cout<<"data c "<<ch<<" bin "<<i<<": "<<dynamic_cast<RooRealVar*>(v_pdfs_roodataset[ch]->get(i)->first())->getVal()<<" weight "<<v_pdfs_roodataset[ch]->weight()<<endl;
		    }
		    f.WriteTObject(&h);

		    TCanvas canvas("can","can",800, 800);
		    canvas.cd();
		    h.Draw("histe");
		    double normS = 0, normB=0;
		    for(int i=0; i<vv_pdfs_norm_varied[ch].size(); i++) {
			    //if(i<v_pdfs_sigproc[ch]) normS+=vv_pdfs_norm_varied[ch][i];
			    //else normB+=vv_pdfs_norm_varied[ch][i];
			    if(i<v_pdfs_sigproc[ch]) normS+=_w_varied->var(vv_pdfs_normNAME[ch][i])->getVal();
			    else  normB+=_w_varied->var(vv_pdfs_normNAME[ch][i])->getVal();
		    }
		    h.SetLineWidth(2);
		    hs.SetLineWidth(2);
		    hb.SetLineWidth(2);
		    hsb.SetLineWidth(2);
		    cout<<"normS = "<<normS<<", normB = "<<normB<<endl;
		    hs.Scale(normS*binwidth);
		    hs.SetLineColor(kBlue);
		    hs.Draw("same");
		    hb.Scale(normB*binwidth);
		    hb.SetLineColor(kRed);
		    hb.Draw("same");
		    hsb.Scale((normS+normB)*binwidth);
		    hsb.SetLineColor(6);
		    hsb.Draw("same");
		    f.WriteTObject(&canvas);

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

    void CountingModel::AddObservedDataSet(string channelname, RooAbsData* rds){
	    int ch =0;
	    for(int c = 0; c<v_pdfs_channelname.size(); c++){
		    if(v_pdfs_channelname[c]==channelname) ch=c;
	    }
	    AddObservedDataSet(ch, rds);
    }
    double CountingModel::EvaluateChi2(double *par, vector< vector< vector<float> > > & vvv_cachPdfValues, int bUseBestEstimateToCalcQ){ 
	    double ret=0;

	    FluctuatedNumbers(par, true, bUseBestEstimateToCalcQ, false);

	    RooRealVar * tmprrv;

	    if(vvv_cachPdfValues.size()==0){
		    vvv_cachPdfValues.resize(vv_pdfs.size());
		    for(int ch=0; ch<vv_pdfs.size(); ch++){
			    vvv_cachPdfValues[ch].resize(vv_pdfs[ch].size());
			    int ntot = int(v_pdfs_roodataset_tmp[ch]->numEntries());
			    for(int p=0; p<vv_pdfs[ch].size();p++){
				    vvv_cachPdfValues[ch][p].resize(ntot);
			    }
		    }
	    }

	    //int maxsets_forcaching = 20;
	    if(vvvv_pdfs_ChProcSetEvtVals.size()==0){
		    vvvv_pdfs_ChProcSetEvtVals.resize(vv_pdfs.size());
		    for(int ch=0; ch<vv_pdfs.size(); ch++){
			    vvvv_pdfs_ChProcSetEvtVals[ch].resize(vv_pdfs[ch].size());
			    for(int p=0; p<vv_pdfs[ch].size();p++){
				    vvvv_pdfs_ChProcSetEvtVals[ch][p].resize(maxsets_forcaching); 
				    int ntot = int(v_pdfs_roodataset_tmp[ch]->numEntries());
				    for(int i=0; i<maxsets_forcaching; i++){
					    vvvv_pdfs_ChProcSetEvtVals[ch][p][i].resize(ntot); 
				    }
			    }
		    }
	    }
	    if(_debug>=100)cout<<" hi0 "<<endl;
	    if(vvvv_pdfs_ChProcSetParVals.size()==0){
		    vvvv_pdfs_ChProcSetParVals.resize(vv_pdfs.size());
		    for(int ch=0; ch<vv_pdfs.size(); ch++){
			    vvvv_pdfs_ChProcSetParVals[ch].resize(vv_pdfs[ch].size());
			    for(int p=0; p<vv_pdfs[ch].size();p++){
				    vvvv_pdfs_ChProcSetParVals[ch][p].resize(maxsets_forcaching); 
				    for(int i=0; i<maxsets_forcaching; i++){
					    vvvv_pdfs_ChProcSetParVals[ch][p][i].resize(vvv_pdfs_nuisancesindex[ch][p].size()); 
					    if(vvv_pdfs_nuisancesindex[ch][p].size()>0)vvvv_pdfs_ChProcSetParVals[ch][p][i][0] = -99998.;
				    }
			    }
		    }
	    }
	    if(_debug>=100)cout<<" hi1 "<<endl;
	    if(vv_pdfs_curSetIndex.size()==0){
		    vv_pdfs_curSetIndex.resize(vv_pdfs.size());
		    for(int ch=0; ch<vv_pdfs.size(); ch++){
			    vv_pdfs_curSetIndex[ch].resize(vv_pdfs[ch].size());
			    for(int p=0; p<vv_pdfs[ch].size(); p++) vv_pdfs_curSetIndex[ch][p]=-1;
		    }
	    }
	    if(_debug>=100)cout<<" hi2 "<<endl;
	    if(TMP_vvpdfs_chprocINT.size()==0){
		    TMP_vvpdfs_chprocINT.resize(vv_pdfs.size());
		    for(int ch=0; ch<vv_pdfs.size(); ch++) {
			    TMP_vvpdfs_chprocINT[ch].resize(vv_pdfs[ch].size());
			    for(int p=0; p<vv_pdfs[ch].size(); p++) TMP_vvpdfs_chprocINT[ch][p]=-1;
		    }
	    }
	    if(maxsets_forcaching>0){
		    for(int ch=0; ch<vv_pdfs.size(); ch++) {
			    for(int p=0; p<vv_pdfs[ch].size(); p++) TMP_vvpdfs_chprocINT[ch][p]=-1;
			    //	for(int p=0; p<vv_pdfs[ch].size(); p++) vv_pdfs_curSetIndex[ch][p]=-1;
		    }
	    }
	    double btot = 0, stot=0, sbtot=0;
	    double tmp=0, tmp2=0, tmp3=0, retch=0;
	    int ntot;
	    double weight=1;
	    int nsigproc = 1;
	    double firstSigProcPdfVal = -1;
	    for(int ch=0; ch<vv_pdfs.size(); ch++){
		    btot = 0; stot=0;
		    tmp=0; retch=0;
		    ntot = int(v_pdfs_roodataset_tmp[ch]->numEntries());
		    for(int i=0; i<vv_pdfs_norm_varied[ch].size(); i++){
			    if(i>=v_pdfs_sigproc[ch]) btot+=vv_pdfs_norm_varied[ch][i]; // FIXME HGG // already multiplied by extra norm in FluctuatedNumbers
			    else stot+=vv_pdfs_norm_varied[ch][i];
		    }
		    sbtot=stot+btot; if(sbtot<=0) {continue;} //cout<<"ERROR: evaluateCh2:  s+b <= 0: "<<sbtot<<endl; exit(1); 
		    //RooArgSet vars(*(_w_varied->var(v_pdfs_observables[ch])));
		    RooArgSet vars(*(v_pdfs_roodataset_real[ch])->get());

		    /*
		       for(int i=0; i<ntot; i++){
		       _w_varied->var(v_pdfs_observables[ch])->setVal(( dynamic_cast<RooRealVar*>(v_pdfs_roodataset_tmp[ch]->get(i)->first()))->getVal());
		       tmp = 0;  ///////////////
		       if(stot!=0) tmp += stot*_w_varied->pdf(v_pdfs_s[ch])->getVal(&vars);  //give some warning message when r=0
		       tmp += btot*_w_varied->pdf(v_pdfs_b[ch])->getVal(&vars);
		       retch -= (tmp>0?log(tmp):0);
		       }
		     */

		    nsigproc = v_pdfs_sigproc[ch];
		    //TString stemp = "";
		    //bool bupdated = false;
		    for(int p=0; p<vv_pdfs[ch].size(); p++) TMP_vvpdfs_chprocINT[ch][p]=-1;
		    for(int i=0; i<ntot; i++){
			    //stemp+=" evt "; stemp+=i;
			    //tmprrv=_w_varied->var(v_pdfs_observables[ch]);
			    //tmprrv->setDirtyInhibit(1);
			    //tmprrv->setVal(( dynamic_cast<RooRealVar*>(v_pdfs_roodataset_tmp[ch]->get(i)->first()))->getVal());

			    int itmpp=0;

			    std::auto_ptr<TIterator> iter(v_pdfs_roodataset_tmp[ch]->get(i)->createIterator());
			    for(RooRealVar * obs= (RooRealVar*)iter->Next(); obs!=NULL; obs=(RooRealVar*)iter->Next()){

				    if(TString(obs->GetName()).BeginsWith("wgttmp_")) continue;
				    _w_varied->var(obs->GetName())->setVal(obs->getVal());
				    //cout<<" Bin "<<i<<": "<<obs->getVal()<<endl;
			    }

			    weight = v_pdfs_roodataset_tmp[ch]->weight();
			    if(_debug>=10 && i==0) cout<<"channel "<<ch<<",  event "<<i<<": weight = "<<weight<<endl;
			    tmp = 0;  
			    firstSigProcPdfVal = -1;
			    for(int p=0; p<vv_pdfs[ch].size();p++){
				    //stemp+=" proc "; stemp+=p;
				    if(vv_pdfs_norm_varied[ch][p]!=0){
					    if(_debug>=100&&i==0){
						    cout<<"DELETEME 1.1: "<<vv_pdfs_norm_varied[ch][p]<<endl;
						    cout<<vvv_cachPdfValues.size()<<endl;
						    cout<<" ch="<<ch<<" p="<<p<<endl;
						    cout<<" Updated ? "<<(vv_pdfs_statusUpdated[ch][p])<<endl;
						    cout<<"v[ch].size="<<vvv_cachPdfValues[ch].size()<<endl;
						    cout<<"v[ch][p].size="<<vvv_cachPdfValues[ch][p].size()<<endl;
					    }
					    if(vv_pdfs_statusUpdated[ch][p]){
						    if( p<nsigproc && firstSigProcPdfVal >= 0 && b_MultiSigProcShareSamePDF){
							    tmp3 = firstSigProcPdfVal;
						    }else{


							    if(i==0 && maxsets_forcaching>0){
								    for(int ti=0; ti<maxsets_forcaching; ti++){
									    for(int tj=0; tj < vvv_pdfs_nuisancesindex[ch][p].size(); tj++){
										    if((par[vvv_pdfs_nuisancesindex[ch][p][tj]])!=vvvv_pdfs_ChProcSetParVals[ch][p][ti][tj]) break;
										    if(tj==vvv_pdfs_nuisancesindex[ch][p].size()-1)TMP_vvpdfs_chprocINT[ch][p] = ti ;
									    }
								    }	
							    }

							    if(TMP_vvpdfs_chprocINT[ch][p]>=0){
								    tmp3 = vvvv_pdfs_ChProcSetEvtVals[ch][p][TMP_vvpdfs_chprocINT[ch][p]][i];
								    //  float tmp4 =	_w_varied->pdf(vv_pdfs[ch][p].c_str())->getVal(&vars);  //give some warning message when r=0
								    //  float tmp5 =	_w_varied->pdf(vv_pdfs[ch][p].c_str())->getVal(&vars);  //give some warning message when r=0
								    if(0){
									    //	cout<<" tmp3="<<tmp3<<" tmp4="<<tmp4<<" tmp4p="<<tmp5<<endl;
									    cout<<" tmp3_0="<<vvvv_pdfs_ChProcSetEvtVals[ch][p][TMP_vvpdfs_chprocINT[ch][p]][0]<<endl;	
									    for(int tii=0; tii<maxsets_forcaching; tii++){
										    cout<<" tii="<<tii<<" tmp3="<<vvvv_pdfs_ChProcSetEvtVals[ch][p][tii][0]<<"  "<<vvvv_pdfs_ChProcSetEvtVals[ch][p][tii][1]<<endl;
									    }

									    int tmpii = TMP_vvpdfs_chprocINT[ch][p];
									    for(int ppi=0; ppi<vvvv_pdfs_ChProcSetParVals[ch][p][tmpii].size(); ppi++)cout<<" TTTT "<< vvvv_pdfs_ChProcSetParVals[ch][p][tmpii][ppi]/*<<"  name="<<v_uncname[vvv_pdfs_nuisancesindex[ch][p][ppi]-1]*/<<" index= "<<vvv_pdfs_nuisancesindex[ch][p][ppi]<<endl;
									    cout<<" PPPPP "<<par[1]<<" name="<<v_uncname[0] <<endl;
									    cout<<" PPPPP "<<par[2]<<" name="<<v_uncname[1] <<endl;
									    cout<<" PPPPP "<<par[3]<<" name="<<v_uncname[2] <<endl;
									    cout<<" PPPPP "<<par[4]<<" name="<<v_uncname[3] <<endl;





									    exit(1); //

								    }
							    }else{	
								    if(i==0 && maxsets_forcaching>0){ // move to next set, and if >= maxsets_forcaching, take the mode 
									    vv_pdfs_curSetIndex[ch][p]+=1; 
									    //cout<<" curSetIndex = "<<vv_pdfs_curSetIndex[ch][p]<<endl;
									    if(vv_pdfs_curSetIndex[ch][p]>=maxsets_forcaching) vv_pdfs_curSetIndex[ch][p]-=maxsets_forcaching;
									    for(int tj=0; tj<vvv_pdfs_nuisancesindex[ch][p].size(); tj++){
										    vvvv_pdfs_ChProcSetParVals[ch][p][vv_pdfs_curSetIndex[ch][p]][tj]=(par[vvv_pdfs_nuisancesindex[ch][p][tj]]);
										    //// make it to not work ,   ---> to test timing
										    //if(tj==vvv_pdfs_nuisancesindex[ch][p].size()-1) vvvv_pdfs_ChProcSetParVals[ch][p][vv_pdfs_curSetIndex[ch][p]][tj]=-9999998.;
									    }
								    }
								    tmp3 =	_w_varied->pdf(vv_pdfs[ch][p].c_str())->getVal(&vars);  //give some warning message when r=0
								    //  cout<<" i "<<i<<" val="<<tmp3<<"   *********ch"<<ch<<"p"<<p<<"  before val="<<vvvv_pdfs_ChProcSetEvtVals[ch][p][vv_pdfs_curSetIndex[ch][p]][i]<<endl;
								    if(maxsets_forcaching>0)vvvv_pdfs_ChProcSetEvtVals[ch][p][vv_pdfs_curSetIndex[ch][p]][i] = tmp3; // update the value of the chosen set 
								    _countPdfEvaluation ++ ;
								    //   cout<<" i "<<i<<" val="<<tmp3<<"   *********ch"<<ch<<"p"<<p<<"  curSetIndex="<<vv_pdfs_curSetIndex[ch][p]<<"  after val="<<vvvv_pdfs_ChProcSetEvtVals[ch][p][vv_pdfs_curSetIndex[ch][p]][i]<<endl;
							    }
							    if(p<nsigproc) firstSigProcPdfVal = tmp3;
						    }
						    //bupdated = true;
						    tmp2 = vv_pdfs_norm_varied[ch][p]*tmp3;
						    vvv_cachPdfValues[ch][p][i]=tmp3;
						    if(_debug>=100&&i==0)cout<<" new: "<<tmp3<<endl;
						    //stemp+=" val "; stemp+=tmp3;
					    }else {
						    if(_debug==102)tmp3 =	_w_varied->pdf(vv_pdfs[ch][p].c_str())->getVal(&vars);  //give some warning message when r=0
						    tmp2=vv_pdfs_norm_varied[ch][p]*vvv_cachPdfValues[ch][p][i];

						    if(_debug==102 && i==0)cout<<" caching comparison: "<<tmp3<<" "<<vvv_cachPdfValues[ch][p][i]<<endl;
					    }
					    tmp +=tmp2;
				    }
				    if(isnan(tmp) || isinf(tmp)) {
					    cout<<" DELETEME * pdf  "<<tmp<<endl;
					    FlagAllChannels();
					    return 10e9;
				    }
			    }//loop over process
			    retch -= (tmp>0?log(tmp):0)*weight;
		    }//loop over events
		    //if(bupdated) cout<<stemp<<endl;

		    retch+=stot;
		    retch+=btot;
		    retch-= (v_pdfs_roodataset_tmp[ch]->sumEntries());
		    if(_debug >= 10) cout<<" DELETEME sumEntries ="<<v_pdfs_roodataset_tmp[ch]->sumEntries()<<endl;;
		    if(_debug>=10){
			    cout<<"EvaluateChi2 in channel ["<<v_pdfs_channelname[ch]<<"]: lnQ= "<<retch<<endl;
			    cout<<"\n model_sb"<<endl;
			    _w_varied->pdf(v_pdfs_sb[ch])->getParameters(*(_w->var(v_pdfs_observables[ch])))->Print("V");
		    }
		    ret+=retch;
	    }// loop over channels
	   if(isnan(ret) || isinf(ret)) cout<<" DELETEME **** pdfs ret = "<<ret<<endl;
	    return ret;
    }

    double CountingModel::EvaluateGL(int ch, double xr){ // deprecated
	    cout<<"Deprecated funtion:  EvaluateGL()"<<endl; exit(1);
	    double ret=0;

	    double btot = 0, stot=0;
	    int ntot = int(v_pdfs_roodataset[ch]->sumEntries());
	    double tmp;
	    for(int i=0; i<vv_pdfs_norm_scaled[ch].size(); i++){ //FIXME HGG
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
		    _w_varied->pdf(v_pdfs_sb[ch])->getParameters(*_w->var(v_pdfs_observables[ch]))->Print("V");
	    }
	    return ret;
    }
    void CountingModel::SetDataForUnbinned(TString filename, bool bRealData){
	    vector<RooAbsData*> data; data.clear();
	    for(int i=0; i<v_pdfs_channelname.size(); i++){
		    data.push_back((RooAbsData*)GetTObject(filename.Data(), v_pdfs_channelname[i]));
	    }
	    SetDataForUnbinned(data, bRealData);
    }
    void CountingModel::SetData(TString filename, TString datahistname, bool bRealData){
	    TH1D * h = (TH1D*)(GetTObject(filename.Data(), datahistname.Data()));
	    if(h->GetNbinsX()!=v_channelname.size()){cout<<"ERROR: file "<<filename<<" hist "<<datahistname<<" only has "<<h->GetNbinsX()
		    <<" != nchannels="<<v_channelname.size()<<endl; exit(1);}
	    vector<double> data; data.clear();
	    for(int i=0; i<v_data.size(); i++){
		    data.push_back(h->GetBinContent(i+1));
	    }
	    SetData(data, bRealData);
    }
    void CountingModel::SetDataForUnbinned(vector< RooAbsData* > data, bool bRealData){
	    /*
	       for(int ch=0; ch<vv_pdfs.size(); ch++){
	       if(v_pdfs_roodataset[ch])delete v_pdfs_roodataset[ch];
	       if(bRealData){
	       if(v_pdfs_roodataset_real[ch])delete v_pdfs_roodataset[ch];
	       }
	       }
	       */
	    v_pdfs_roodataset.clear();
	    v_pdfs_roodataset_tmp.clear();
	    if(bRealData){
		    v_pdfs_roodataset_real.clear();
	    }
	    for(int ch=0; ch<vv_pdfs.size(); ch++){
		    vector<double> vtmp;vtmp.clear();
		    v_pdfs_roodataset.push_back(data[ch]);
		    v_pdfs_roodataset_tmp.push_back(data[ch]);
		    if(bRealData) v_pdfs_roodataset_real.push_back(data[ch]);

		    /*
		       for(int i=0; i<v_pdfs_roodataset[ch]->sumEntries(); i++){
		       RooRealVar *r = dynamic_cast<RooRealVar*>(v_pdfs_roodataset[ch]->get(i)->first());
		       vtmp.push_back( r->getVal() );
		       }
		       vv_pdfs_data.push_back(vtmp);
		       */
	    }
    }
    void CountingModel::SetTmpDataForUnbinned(vector< RooAbsData* > data){
	    v_pdfs_roodataset_tmp.clear();
	    for(int ch=0; ch<vv_pdfs.size(); ch++){
		    v_pdfs_roodataset_tmp.push_back(data[ch]);
	    }
    }
    void CountingModel::SetToyForUnbinned(vector< RooAbsData* > data){
	    v_pdfs_roodataset_toy.clear();
	    for(int ch=0; ch<vv_pdfs.size(); ch++){
		    v_pdfs_roodataset_toy.push_back(data[ch]);
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
	    if(b_ForceSymmetryError) {
		    if (uncertainty_in_relative_fraction_up > uncertainty_in_relative_fraction_down) uncertainty_in_relative_fraction_down=uncertainty_in_relative_fraction_up;
		    if(uncertainty_in_relative_fraction_down > uncertainty_in_relative_fraction_up ) uncertainty_in_relative_fraction_up = uncertainty_in_relative_fraction_down;
	    }
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
		    v_uncFloatInFit.push_back(1);
	    }
	    if(_debug>=10) cout<<" AddUncertaintyOnShapeNorm for "<<index_channel<<"th channel, "<<index_sample<<"th proc,  unc ["<<uncname<<"]"<<endl;
	    AddUncertaintyOnShapeNorm(index_channel, index_sample, uncertainty_in_relative_fraction_down, uncertainty_in_relative_fraction_up, pdf_type, index_correlation );
    }
    void CountingModel::AddUncertaintyOnShapeNorm(string c, int index_sample, double uncertainty_in_relative_fraction_down, double uncertainty_in_relative_fraction_up, int pdf_type, string uncname ){
	    int index_channel = -1;
	    for(int i=0; i<v_pdfs_channelname.size(); i++){
		    if(v_pdfs_channelname[i]==c) index_channel=i;
	    }
	    AddUncertaintyOnShapeNorm(index_channel, index_sample, uncertainty_in_relative_fraction_down, uncertainty_in_relative_fraction_up, pdf_type, uncname);

    }

    void CountingModel::AddUncertaintyOnShapeNorm(int index_channel, int index_sample, double rho, double rho_err, double B, int pdf_type, string uncname ){
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
		    v_uncFloatInFit.push_back(1);
	    }
	    AddUncertaintyOnShapeNorm(index_channel, index_sample, rho, rho_err, B, pdf_type, index_correlation );
    }
    void CountingModel::AddUncertaintyOnShapeNorm(string c, int index_sample, double rho, double rho_err, double B, int pdf_type, string uncname ){
	    int index_channel = -1;
	    for(int i=0; i<v_pdfs_channelname.size(); i++){
		    if(v_pdfs_channelname[i]==c) index_channel=i;
	    }
	    AddUncertaintyOnShapeNorm(index_channel, index_sample, rho, rho_err, B, pdf_type, uncname);
    }

    // From SideBand
    // when B is large, it's close to Gaussian. then use Gaussian, don't use Gamma function
    void CountingModel::AddUncertaintyOnShapeNorm(int index_channel, int index_sample, double rho, double rho_err, double B, int pdf_type, int index_correlation ){
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
	    vvv_pdfs_normvariation.at(index_channel).at(index_sample).push_back(vunc);
	    vvv_pdfs_pdftype.at(index_channel).at(index_sample).push_back(pdf_type);
	    vvv_pdfs_idcorrl.at(index_channel).at(index_sample).push_back(index_correlation);
    }
    double CountingModel::EvaluateGL(vector< vector<double> > vnorms, vector<double> vparams, double xr, vector< vector<double> > & vvs, vector< vector<double> > &vvb){ 
	    double ret=0;

	    if(vvs.size()>0 && vvb.size()>0){
		    for(int ch=0; ch<vv_pdfs.size(); ch++){
			    int ntot = int(v_pdfs_roodataset[ch]->numEntries());
			    double tmp=0;
			    for(int i=0; i<ntot; i++){
				    v_pdfs_roodataset[ch]->get(i);
				    tmp+= log( xr*vvs[ch][i] + vvb[ch][i])*v_pdfs_roodataset[ch]->weight();
			    }
			    ret+=tmp;
		    }
	    }else{
		    if(_debug>=100){cout<<"In EvaluateGL:  xr ="<<xr<<endl;};
		    for(int i=0; i<v_pdfs_floatParamsName.size(); i++){
			    if(_debug>=100) cout<<" EvaluateGL: setting value of parameter ["<<v_pdfs_floatParamsName[i]<<"] = "<<vparams[i]<<endl;
			    _w_varied->var(v_pdfs_floatParamsName[i].c_str())->setVal(vparams[i]);
		    }
		    double tmps, tmpb;
		    vector<double> vs, vb;
		    for(int ch=0; ch<vv_pdfs.size(); ch++){
			    vs.clear(); vb.clear();
			    double btot = 0, stot=0;
			    int ntot = int(v_pdfs_roodataset[ch]->numEntries());
			    double tmp=0;
			    for(int i=0; i<vnorms[ch].size(); i++){
				    if(i>=v_pdfs_sigproc[ch]) btot+=vnorms[ch][i];
				    else stot+=vnorms[ch][i];

				    _w_varied->var(vv_pdfs_normNAME[ch][i])->setVal(vnorms[ch][i]);
				    if(_debug>=100)cout<<vv_pdfs_normNAME[ch][i]<<" "<<vnorms[ch][i]<<endl;
			    }
			    //RooArgSet vars(*(_w->var(v_pdfs_observables[ch])));
			    RooArgSet vars(*(v_pdfs_roodataset_real[ch])->get());
			    for(int i=0; i<ntot; i++){
				    //_w_varied->var(v_pdfs_observables[ch])->setVal(( dynamic_cast<RooRealVar*>(v_pdfs_roodataset[ch]->get(i)->first()))->getVal());

				    std::auto_ptr<TIterator> iter(v_pdfs_roodataset[ch]->get(i)->createIterator());
				    for(RooRealVar * obs= (RooRealVar*)iter->Next(); obs!=0; obs=(RooRealVar*)iter->Next()){
					    if(TString(obs->GetName()).BeginsWith("wgttmp_")) continue;
					    _w_varied->var(obs->GetName())->setVal(obs->getVal());
				    }
											              
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
				    tmps =  stot*_w_varied->pdf(v_pdfs_s[ch])->getVal(&vars); tmpb =  btot*_w_varied->pdf(v_pdfs_b[ch])->getVal(&vars);
				    vs.push_back(tmps); vb.push_back(tmpb);
				    tmp+= log(xr*tmps+tmpb) * v_pdfs_roodataset[ch]->weight();
				    //tmp+= log( xr*stot*_w_varied->pdf(v_pdfs_s[ch])->getVal(&vars)
				    //		+btot*_w_varied->pdf(v_pdfs_b[ch])->getVal(&vars));
			    }

			    if(_debug>=100){
				    cout<<"EvaluateGL in channel ["<<v_pdfs_channelname[ch]<<"]: gl = "<<tmp<<endl;
				    cout<<"\n model_sb"<<endl;
				    _w_varied->pdf(v_pdfs_sb[ch])->getParameters(*_w->var(v_pdfs_observables[ch]))->Print("V");
			    }
			    ret+=tmp;
			    vvs.push_back(vs); vvb.push_back(vb);
		    }
	    }
	    if(_debug>=100) {
		    cout<<" Unbinned Model EvaluateGL all channels = " << ret <<endl;
	    }
	    return ret;
    }
   bool CountingModel::AddUncertaintyOnShapeParam(string pname){ // add flatParam   // taken from workspace
	    if(_debug>=10)cout<<"In AddUncertaintyOnShapeParam  Adding floating parameter with flatParam: "<<pname<<endl;
	    if(_w_varied->var(pname.c_str())==NULL) {
		    cout<<" parameter "<<pname<<" not exist in the added channels,   skip it"<<endl;
		    if(_w_varied->arg(pname.c_str())) cout<<" but it exist as RooAbsArg"<<endl;
		    _w_varied->Print();
		    return false;
	    }

	    int index_correlation = -1; // numeration starts from 1
	    for(int i=0; i<v_uncname.size(); i++){
		    if(v_uncname[i]==pname){
			    index_correlation = i+1;	
			    cout<<" There are two shape parameters with same name = "<<pname<<", skip it"<<endl;
			    return false;
		    }
	    }
	    if(index_correlation<0)  {
		    index_correlation = v_uncname.size()+1;
		    v_uncname.push_back(pname);
		    v_uncFloatInFit.push_back(1);
	    }

	    RooRealVar* rrv = (RooRealVar*)_w_varied->var(pname.c_str());
	    double norminalValue = rrv->getVal();
	    double rangeMin = rrv->getMin();
	    double rangeMax = rrv->getMax();
	    if(_debug) {
		    rrv->Print();
		    cout<<norminalValue<<", ["<<rangeMin<<", "<<rangeMax<<"]"<<endl;
	    }


	    vector<double> vunc; vunc.clear(); 
	    vunc.push_back(norminalValue);
	    vunc.push_back(rangeMax-rangeMin);
	    vunc.push_back(0);
	    vunc.push_back(rangeMin);
	    vunc.push_back(rangeMax);
	    if(v_pdfs_floatParamsUnc.size() <= index_correlation){
		    for(int i = v_pdfs_floatParamsUnc.size(); i<=index_correlation; i++){
			    v_pdfs_floatParamsUnc.push_back(vunc);
		    }
	    }else v_pdfs_floatParamsUnc[index_correlation] = vunc;
	    if(_debug) cout<<" v_pdfs_floatParamsUnc.size() = "<<v_pdfs_floatParamsUnc.size()<<" index_correlation="<<index_correlation<<endl;
	    if(_debug)cout<<" * Adding flat floating parameter: "<<pname<<endl;

	    TString s = pname;
	    v_pdfs_floatParamsName.push_back(pname);
	    v_pdfs_floatParamsIndcorr.push_back(index_correlation);
	    v_pdfs_floatParamsType.push_back(typeFlat);

	    if(_debug)cout<<"FlatParam "<<pname<<": nominal value = "<<norminalValue<<" ["<<rangeMin<<","<<rangeMax<<"]"<<endl;

	    return true;
    }
    bool CountingModel::AddFlatParam(string pname, double norminalValue, double rangeMin, double rangeMax){ // add flatParam,  taken from text file , ie. data cards
	    if(_debug>=10)cout<<"In AddUncertaintyOnShapeParam  Adding floating parameter with flatParam: "<<pname<<endl;
	    if(_w_varied->var(pname.c_str())==NULL) {
		    cout<<" parameter "<<pname<<" not exist in the added channels,   skip it"<<endl;
		    if(_w_varied->arg(pname.c_str())) cout<<" but it exist as RooAbsArg"<<endl;
		    _w_varied->Print();
		    return false;
	    }

	    int index_correlation = -1; // numeration starts from 1
	    for(int i=0; i<v_uncname.size(); i++){
		    if(v_uncname[i]==pname){
			    index_correlation = i+1;	
			    cout<<" There are two shape parameters with same name = "<<pname<<", skip it"<<endl;
			    return false;
		    }
	    }
	    if(index_correlation<0)  {
		    index_correlation = v_uncname.size()+1;
		    v_uncname.push_back(pname);
		    v_uncFloatInFit.push_back(1);
	    }

	    RooRealVar* rrv = (RooRealVar*)_w_varied->var(pname.c_str());
/*
	    double norminalValue = rrv->getVal();
	    double rangeMin = rrv->getMin();
	    double rangeMax = rrv->getMax();
*/	    if(_debug) {
		    rrv->Print();
		    cout<<norminalValue<<", ["<<rangeMin<<", "<<rangeMax<<"]"<<endl;
	    }


	    vector<double> vunc; vunc.clear(); 
	    vunc.push_back(norminalValue);
	    vunc.push_back(rangeMax-rangeMin);
	    vunc.push_back(0);
	    vunc.push_back(rangeMin);
	    vunc.push_back(rangeMax);
	    if(v_pdfs_floatParamsUnc.size() <= index_correlation){
		    for(int i = v_pdfs_floatParamsUnc.size(); i<=index_correlation; i++){
			    v_pdfs_floatParamsUnc.push_back(vunc);
		    }
	    }else v_pdfs_floatParamsUnc[index_correlation] = vunc;
	    if(_debug) cout<<" v_pdfs_floatParamsUnc.size() = "<<v_pdfs_floatParamsUnc.size()<<" index_correlation="<<index_correlation<<endl;
	    if(_debug)cout<<" * Adding flat floating parameter: "<<pname<<endl;

	    TString s = pname;
	    v_pdfs_floatParamsName.push_back(pname);
	    v_pdfs_floatParamsIndcorr.push_back(index_correlation);
	    v_pdfs_floatParamsType.push_back(typeFlat);

	    if(_debug)cout<<"FlatParam "<<pname<<": nominal value = "<<norminalValue<<" ["<<rangeMin<<","<<rangeMax<<"]"<<endl;

	    return true;
    }
    bool CountingModel::AddUncertaintyOnShapeParam(string pname, double mean, double sigmaL, double sigmaR, double rangeMin, double rangeMax ){
	    if(_debug>=10)cout<<"In AddUncertaintyOnShapeParam  Adding floating parameter: "<<pname<<endl;
	    if(_w_varied->var(pname.c_str())==NULL) {
		    cout<<" parameter "<<pname<<" not exist in the added channels,   skip it"<<endl;
		    return false;
	    }
	    if(_debug) {
		    RooRealVar * rrv = _w_varied->var(pname.c_str());
		    cout<<" parameter "<<pname<<" in workspace: "<<endl;
		    rrv->Print("");
		    if(_debug>=10)rrv->Print("V");
		    cout<<"\n"<<endl;
		    cout<<" in datacard: "<<endl;
		    cout<<" mean = "<<mean<<endl;
		    cout<<" sigma L= "<<sigmaL<<", R="<<sigmaR<<endl;
		    cout<<" range = ["<<rangeMin<<","<<rangeMax<<"]"<<endl;
	    }

	    int index_correlation = -1; // numeration starts from 1
	    sigmaL = fabs(sigmaL);
	    sigmaR = fabs(sigmaR);
	    for(int i=0; i<v_uncname.size(); i++){
		    if(v_uncname[i]==pname){
			    index_correlation = i+1;	
			    cout<<" There are two shape parameters with same name = "<<pname<<", skip it"<<endl;
			    return false;
		    }
	    }
	    if(index_correlation<0)  {
		    index_correlation = v_uncname.size()+1;
		    v_uncname.push_back(pname);
		    v_uncFloatInFit.push_back(1);
	    }

	    if(rangeMax==rangeMin){
		    rangeMin = mean - 4*sigmaL;
		    rangeMax = mean + 4*sigmaR;
	    }

	    RooRealVar * rrv = _w_varied->var(pname.c_str());
	    rrv->setRange(rangeMin, rangeMax);

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
	    if(_debug) cout<<" v_pdfs_floatParamsUnc.size() = "<<v_pdfs_floatParamsUnc.size()<<" index_correlation="<<index_correlation<<endl;

	    TString s = pname;
	    if(_debug)cout<<" * Adding bifurcated floating parameter: "<<pname<<endl;
	    _w_varied->factory(TString::Format("%s_x[%f,%f]", pname.c_str(), rangeMin, rangeMax));
	    _w_varied->factory(TString::Format("%s_mean[%f,%f,%f]", pname.c_str(), mean, rangeMin, rangeMax));
	    if(_debug)cout<<pname<<"_x added"<<endl;
	    _w_varied->factory(TString::Format("BifurGauss::%s_bfg(%s_x, %s_mean, %f, %f )", pname.c_str(), pname.c_str(), pname.c_str(), sigmaL, sigmaR));
	    if(_debug)cout<<pname<<"_bfg added"<<endl;

	    v_pdfs_floatParamsName.push_back(pname);
	    v_pdfs_floatParamsIndcorr.push_back(index_correlation);
	    v_pdfs_floatParamsType.push_back(typeBifurcatedGaussian);

	    return true;
    }

    void CountingModel::AddUncertaintyAffectingShapeParam(string uname, string pname, double sigmaL, double sigmaR ){

	    // for the uncertainty uname,  assign it a range for truncation ....     sigmaL sigmaR are absolute values, need to translate to  gaussian unit
	    //
	    // type should be gaussain ,  truncated
	    //
	    // for each uncertainty source:  need a range,  and a list of parameters which it affects 

	    bool pname_added = false;
	    for(int i=0; i<v_uncname.size(); i++){
		    if(v_uncname[i]==pname){
			    pname_added = true;
		    }
	    }
	    if(!pname_added){
		    cout<<"In AddUncertaintyAffectingShapeParam: parameter["<<pname<<"] not added yet"<<endl;
		    cout<<"Need a line in data card: \""<<pname<<"\" param mean sigma"<<endl;
		    exit(1);
	    }


	    int index_correlation = -1; // numeration starts from 1
	    sigmaL = fabs(sigmaL);
	    sigmaR = fabs(sigmaR);
	    for(int i=0; i<v_uncname.size(); i++){
		    if(v_uncname[i]==uname){
			    index_correlation = i+1;	
			    break;
		    }
	    }
	    if(index_correlation<0)  {
		    index_correlation = v_uncname.size()+1;
		    v_uncname.push_back(uname);
		    v_uncFloatInFit.push_back(1);
	    }


	    vector<double> v; v.push_back(float(index_correlation)); v.push_back(sigmaL); v.push_back(sigmaR);
	    MapStrVV::iterator iter = map_param_sources.find(pname);
	    if (iter != map_param_sources.end() ) {
		    (iter->second).push_back(v);
	    }else{
		    vector< vector<double> > vv; vv.push_back(v);
		    map_param_sources[pname]=vv;
	    }
	    if(_debug)cout<<"DELETEME: map_param_sources.size = "<<map_param_sources.size()<<endl;
    }
    const vector< RooAbsData* >& CountingModel::Get_v_pdfs_roodataset(){
	    //cout<<"in Get_v_pdfs_roodataset"<<endl;
	    //cout<<" v_pdfs_roodataset.size =  "<<v_pdfs_roodataset.size()<<endl;
	    //v_pdfs_roodataset[0]->Print();
	    return v_pdfs_roodataset;
    }; // in each channel, it has a list of events
    void CountingModel::SetMass(double m){
	    if(_w and _w_varied) {
		    if(_w->var("MH")){
			    cout<<" ** running with Mass point = "<<m<<endl;
			    if(m<=0) { cout<<"ERROR: There is variable MH in workspace, but the input mass is "<<m<<endl; exit(1); };
			    _w->var("MH")->setVal(m);
			    _w_varied->var("MH")->setVal(m);
		    }
	    }else{
		    cout<<"ERROR: SetMass() must be invoked after ConfigureModel"<<endl; exit(1);
	    }
    }
    void CountingModel::FlagChannelsWithParamsUpdated(int i){
	    if(_debug>=100 )cout<<" DELETEME FlagChannelsWithParamsUpdated "<<(i>0?v_uncname[i-1]:"r")<<endl;
	    if(vv_pdfs_statusUpdated.size()==0){
		    vv_pdfs_statusUpdated.resize(vv_pdfs.size());
		    for(int ch=0; ch<vv_pdfs.size(); ch++) vv_pdfs_statusUpdated[ch].resize(vv_pdfs[ch].size());
		    for(int j=0; j<vvp_pdfs_connectNuisBinProc[i].size(); j++){
			    vv_pdfs_statusUpdated[vvp_pdfs_connectNuisBinProc[i][j].first][vvp_pdfs_connectNuisBinProc[i][j].second]=true;
		    }
	    }
	    if(vv_statusUpdated.size()==0){
		    vv_statusUpdated.resize(vv_exp_sigbkgs.size());
		    for(int ch=0; ch<vv_exp_sigbkgs.size(); ch++) vv_statusUpdated[ch].resize(vv_exp_sigbkgs[ch].size());
		    for(int j=0; j<vvp_connectNuisBinProc[i].size(); j++){
			    vv_statusUpdated[vvp_connectNuisBinProc[i][j].first][vvp_connectNuisBinProc[i][j].second]=true;
		    }
	    }
	    for(int j=0; j<vvp_pdfs_connectNuisBinProc[i].size(); j++){
		    vv_pdfs_statusUpdated[vvp_pdfs_connectNuisBinProc[i][j].first][vvp_pdfs_connectNuisBinProc[i][j].second]=true;
	    }
	    for(int j=0; j<vvp_connectNuisBinProc[i].size(); j++){
		    if(_debug>100){
			    cout<<"DELETEME "<<i << " "<<j<<endl;
			    cout<<"ch "<<vvp_connectNuisBinProc[i][j].first<<" proc "<<vvp_connectNuisBinProc[i][j].second<<endl;;
			    cout<<" number of chs "<<vv_statusUpdated.size()<<endl;
		    }
		    vv_statusUpdated[vvp_connectNuisBinProc[i][j].first][vvp_connectNuisBinProc[i][j].second]=true;
	    }
	    if(_debug>=100 )cout<<" DELETEME end FlagChannelsWithParamsUpdated "<<(i>0?v_uncname[i-1]:"r")<<endl;
    }
    void CountingModel::UnFlagAllChannels(bool b){
	    if(_debug>=100) cout<<"in UnFlagAllChannels b="<<b<<endl;
	    if(1){
		    for(int ch=0; ch<vv_pdfs_statusUpdated.size(); ch++){
			    for(int p=0; p<vv_pdfs_statusUpdated[ch].size(); p++)
				    vv_pdfs_statusUpdated[ch][p]=false;
		    }
		    for(int ch=0; ch<vv_statusUpdated.size(); ch++){
			    for(int p=0; p<vv_statusUpdated[ch].size(); p++)
				    vv_statusUpdated[ch][p]=false;
		    }
	    }else{
		    for(int i=0; i<vvp_connectNuisBinProc.size(); i++){
			    for(int j=0; j<vvp_pdfs_connectNuisBinProc[i].size(); j++){
				    vv_pdfs_statusUpdated[vvp_pdfs_connectNuisBinProc[i][j].first][vvp_pdfs_connectNuisBinProc[i][j].second]=false;
			    }
			    for(int j=0; j<vvp_connectNuisBinProc[i].size(); j++){
				    cout<<"DELETEME "<<i << " "<<j<<endl;
				    cout<<"ch "<<vvp_connectNuisBinProc[i][j].first<<" proc "<<vvp_connectNuisBinProc[i][j].second<<endl;;
				    vv_statusUpdated[vvp_connectNuisBinProc[i][j].first][vvp_connectNuisBinProc[i][j].second]=false;
			    }
		    }
	    }
	    if(_debug>=100) cout<<"done in UnFlagAllChannels b="<<b<<endl;
    }
    void CountingModel::FlagAllChannels(){
	    if(_debug>=100) cout<<"in FlagAllChannels"<<endl;
	    vv_pdfs_statusUpdated.resize(vv_pdfs.size());
	    for(int ch=0; ch<vv_pdfs.size(); ch++) vv_pdfs_statusUpdated[ch].resize(vv_pdfs[ch].size());
	    for(int ch=0; ch<vv_pdfs_statusUpdated.size(); ch++){
		    for(int p=0; p<vv_pdfs_statusUpdated[ch].size(); p++)
			    vv_pdfs_statusUpdated[ch][p]=true;
	    }
	    vv_statusUpdated.resize(vv_exp_sigbkgs.size());
	    for(int ch=0; ch<vv_exp_sigbkgs.size(); ch++) vv_statusUpdated[ch].resize(vv_exp_sigbkgs[ch].size());
	    for(int ch=0; ch<vv_statusUpdated.size(); ch++){
		    for(int p=0; p<vv_statusUpdated[ch].size(); p++)
			    vv_statusUpdated[ch][p]=true;
	    }
    }
    void CountingModel::SetFlatParameterRange(int id, double middle, double errLow, double errUp){
	    if(id>= v_pdfs_floatParamsUnc.size()) return;
	    for(int ip=0; ip<v_pdfs_floatParamsIndcorr.size(); ip++){
		    if(v_pdfs_floatParamsIndcorr[ip]==id) {
			    if(v_pdfs_floatParamsType[ip]!=typeFlat) return;
		    }
	    }

	    // for flatParam: vector-->    0  norminal value, 1 max-min,  2 dumy, 3 min, 4 max 
	    if(_debug>=10)cout<<" beginning: middle "<<v_pdfs_floatParamsUnc[id][0]<<", delta "<<v_pdfs_floatParamsUnc[id][1]<<", min "<<v_pdfs_floatParamsUnc[id][3]<<", max "<<v_pdfs_floatParamsUnc[id][4]<<endl;
	    v_pdfs_floatParamsUnc[id][0] = middle;
	    if(errLow>v_pdfs_floatParamsUnc[id][3]) v_pdfs_floatParamsUnc[id][3]=errLow;
	    if(errUp<v_pdfs_floatParamsUnc[id][4]) v_pdfs_floatParamsUnc[id][4]=errUp;
	    v_pdfs_floatParamsUnc[id][1] = v_pdfs_floatParamsUnc[id][4]-v_pdfs_floatParamsUnc[id][3];
	    if(_debug>=10)cout<<" ending: middle "<<v_pdfs_floatParamsUnc[id][0]<<", delta "<<v_pdfs_floatParamsUnc[id][1]<<", min "<<v_pdfs_floatParamsUnc[id][3]<<", max "<<v_pdfs_floatParamsUnc[id][4]<<endl;
    }
    void CountingModel::SetFlatNormalizationRange(int id, double errLow, double errUp){
	    for(int c=0; c< vvv_pdftype.size(); c++){
		    for(int s=0; s<vvv_pdftype[c].size(); s++){
			    for(int u=0; u<vvv_pdftype[c][s].size(); u++){
				    if(vvv_idcorrl[c][s][u]==id && vvv_pdftype[c][s][u]!=typeFlat) return;
				    if(vvv_idcorrl[c][s][u]==id && vvv_pdftype[c][s][u]==typeFlat) {
					    double range = vvvv_uncpar[c][s][u][1]-vvvv_uncpar[c][s][u][0];
					    double oldmin = vvvv_uncpar[c][s][u][0];
					    double oldmax = vvvv_uncpar[c][s][u][1];
					    if(errLow>0) vvvv_uncpar[c][s][u][0]=errLow*range + oldmin;
					    if(errUp>0 && (errUp>=errLow) && errUp<1 ) vvvv_uncpar[c][s][u][1]=errUp*range + oldmin;
				    }
			    }
		    }
	    }

    }
    void CountingModel::DumpFitResults(double *pars, TString ssave){
	VChannelVSample vv = FluctuatedNumbers(pars);	
	vector< vector<TGraphAsymmErrors*> > vvTGraph; vvTGraph.clear();
	for(int c=0; c<vv.size(); c++){
		TString cname = v_channelname[c];
		if(cname.BeginsWith("TH1F_")){
			vector<string> vs;
			StringSplit(vs,v_channelname[c], "_");
			if(vs.size()<5) continue;
			double bincenter = TString(vs[vs.size()-2]).Atof();
			double binlow = TString(vs[vs.size()-3]).Atof();
			double binhigh = TString(vs[vs.size()-1]).Atof() + binlow;
			cname = "";
			for(int i=1; i<vs.size()-3; i++){
				cname += vs[i];
			}
			bool newchannel = false;
			if(c==0){
				newchannel = true;
			}else{
				TString cname1 = v_channelname[c-1];
				if(cname1.BeginsWith("TH1F_")){
					vector<string> vs1;
					StringSplit(vs1,v_channelname[c-1], "_");
					if(vs1.size()<5) { newchannel = true; continue; }
					cname1 = "";
					for(int i=1; i<vs1.size()-3; i++){
						cname1 += vs1[i];
					}
					if(cname != cname1) newchannel = true;
				}else{ newchannel = true; }
			}
			if(newchannel){
				vector<TGraphAsymmErrors*> vTGraph; vTGraph.clear();
				for(int p=0; p<vv[c].size(); p++){
					TGraphAsymmErrors* gr = new TGraphAsymmErrors();
					gr->Set(0);
					vTGraph.push_back(gr);
					TString histname = cname; histname+="_"; histname+=vv_procname[c][p];
					gr->SetName(histname);
				}	
				vvTGraph.push_back(vTGraph);
			}
			for(int p=0; p<vv[c].size(); p++){
				TGraphAsymmErrors * gr = vvTGraph.back()[p];
				gr->Set(gr->GetN()+1);
				gr->SetPoint(gr->GetN()-1, bincenter, vv[c][p]);
				gr->SetPointError(gr->GetN()-1, bincenter-binlow, binhigh-bincenter, 0, 0);
			}
		}
	}
	TFile *f = new TFile(ssave+".root", "RECREATE");
	for(int c=0; c<vvTGraph.size(); c++){
		for(int p=0; p<vvTGraph[c].size(); p++){
			vvTGraph[c][p]->Sort();
			f->WriteTObject(vvTGraph[c][p]);
		}
	}
	f->Close();
    }
};
