//a make the PlotUtilities more configurable, i.e. modulize it as a bunch of classes aa

#include <vector>
#include <iostream>
#include <iomanip>
#include <ctime>
#include <math.h>

#include <TH1F.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TError.h>
#include <TFile.h>
#include <TTree.h>
#include <TPaveText.h>

//#include "config.C"
#include "FloridaStyle.C"
#include "PlotUtilities.h"


#include "TTree.h"
#include "TSystem.h"
// FIXME make a few trees which save all important numbers, so we can easily remake plots 
// and parallel running .....,  combine them at the final step

using std::cout;
using std::endl;
using namespace stats;
char ctmp[255];
string ssave;

const int nReviewPoints=20;
double xReviewPoints[nReviewPoints]={
	115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175, 180, 185, 190, 195, 200, 205, 250
};
string schannels ="HZZ, 4#mu + 4e + 2e2#mu";



/*
string s_log_fullysys = "HZZnoteXcheck/combined/count/log200pb/200pbFull_";
string s_log_nosys = "HZZnoteXcheck/combined/count/log200pb/200pbnosys_";
string s_log_fullsys_correctedM2lnQ = "HZZnoteXcheck/combined/count/log200pb/200pbsys_correctM2lnQ_";
TString s_plots_dir = "plots/HZZ10TeV200pb";
string s_log_4e = "HZZnoteXcheck/singlechannel/4e_200pbsys_";
string s_log_4m = "HZZnoteXcheck/singlechannel/4m_200pbsys_";
string s_log_2e2m = "HZZnoteXcheck/singlechannel/2e2m_200pbsys_";
string s_log_4e_nosys = "HZZnoteXcheck/singlechannel/4e_200pbnosys_";
string s_log_4m_nosys = "HZZnoteXcheck/singlechannel/4m_200pbnosys_";
string s_log_2e2m_nosys = "HZZnoteXcheck/singlechannel/2e2m_200pbnosys_";
TString s_text_topleft = "CMS preliminary: projection for 10TeV, 200 pb^{-1}";
TString s_log_signi = "HZZnoteXcheck/combined/count/sig_log_1fb";  // significances are non sense for 200pb
TString s_log_signi_nosys = "HZZnoteXcheck/combined/count/sig_log_1fb_nosys";
TString s_log_shape_nosys = "HZZnoteXcheck/combined/shape/shape_limits_200pb_nosys";
TString s_log_shape_sys = "HZZnoteXcheck/combined/shape/shape_limits_200pb_sys";

 */
   string s_log_fullysys = "HZZnoteXcheck/combined/count/log1fb/log_";
   string s_log_nosys = "HZZnoteXcheck/combined/count/log1fb/nosys_";
   string s_log_fullsys_correctedM2lnQ = "HZZnoteXcheck/combined/count/log1fb/1fbsys_correctM2lnQ_";
   TString s_plots_dir = "plots/HZZ10TeV1fb";
   string s_log_4e = "HZZnoteXcheck/singlechannel/4e_1fbsys_";
   string s_log_4m = "HZZnoteXcheck/singlechannel/4m_1fbsys_";
   string s_log_2e2m = "HZZnoteXcheck/singlechannel/2e2m_1fbsys_";
   string s_log_4e_nosys = "HZZnoteXcheck/singlechannel/4e_1fbnosys_";
   string s_log_4m_nosys = "HZZnoteXcheck/singlechannel/4m_1fbnosys_";
   string s_log_2e2m_nosys = "HZZnoteXcheck/singlechannel/2e2m_1fbnosys_";
   TString s_log_signi = "HZZnoteXcheck/combined/count/sig_log_1fb";
   TString s_log_signi_nosys = "HZZnoteXcheck/combined/count/sig_log_1fb_nosys";
   TString s_text_topleft = "CMS preliminary: projection for 10TeV, 1 fb^{-1}";
   TString s_log_shape_nosys = "HZZnoteXcheck/combined/shape/shape1fb_nosys";
   //TString s_log_shape_sys = "HZZnoteXcheck/combined/shape/shape1fb_sys";
   TString s_log_shape_sys = "HZZnoteXcheck/combined/count/1fbcount_0.3uncOnSig";

bool bProcessMultiChannels = true;
double projectingXmin = 110; 
double projectingXmax = 260;
double projectingCLsYmin = 1.e-2;
double projectingCLsYmax = 1;
bool   projectingCLsLogY = true;
TString projectingCLsSavePath=s_plots_dir+"/CLs";
string projectingCLsXYtitles="; Higgs mass, m_{H} (GeV/c^{2}); CL_{s}";	
double projectingM2lnQYmin= -20;
double projectingM2lnQYmax=  20;
bool   projectingM2lnQLogY= false;
TString projectingM2lnQSavePath=s_plots_dir+"/M2lnQ";
string projectingM2lnQXYtitles="; Higgs mass, m_{H} (GeV/c^{2}); -2lnQ";	
double projectingRLimitYmin = 0.5;
double projectingRLimitYmax = 30;
bool projectingRLimitLogY = true;
TString projectingRLimitSavePath=s_plots_dir+"/Limit";
string projectingRLimitXYtitles="; Higgs mass, m_{H} (GeV/c^{2}); r=#sigma_{95% CL}/#sigma_{SM}";
double projectingSignificanceYmin = 0.;
double projectingSignificanceYmax = 6;
bool projectingSignificanceLogY = false;
TString projectingSignificanceSavePath=s_plots_dir+"/Significance";
string projectingSignificanceXYtitles="; Higgs mass, m_{H} (GeV/c^{2}); Significance";
TString s_text_topright = "Apr 23 2010";


int debug = 0;
int main(int argc, const char* argv[]){




	sprintf(ctmp, "mkdir %s -p", s_plots_dir.Data() );
	gSystem->Exec(ctmp);

	if(!debug)gErrorIgnoreLevel=5000; // do not show the message when saving plots
	FloridaStyle(); // plots style

	// -----------timing 
	time_t start_time, cur_time;
	time(&start_time);


	vector<double> vr;
	vector<double> vcls;
	vector<double> vp;



	double tmp_cls_r95_nosys[nReviewPoints][5];
	double tmp_cls_r95[nReviewPoints][5];
	double tmp_m2lnQ_b[nReviewPoints][5];
	double tmp_m2lnQ_sb[nReviewPoints][5];
	double tmp_bys_r95[nReviewPoints][6];
	double tmp_bys_r95_nosys[nReviewPoints][6];
	double tmp_cls_cls[nReviewPoints][5];

	double *cls_4e=new double[nReviewPoints];
	double *cls_4m=new double[nReviewPoints];
	double *cls_2e2m=new double[nReviewPoints];

	double *clsr95_4e=new double[nReviewPoints];
	double *clsr95_4m=new double[nReviewPoints];
	double *clsr95_2e2m=new double[nReviewPoints];
	double *bysr95_4e_wos=new double[nReviewPoints];
	double *bysr95_4m_wos=new double[nReviewPoints];
	double *bysr95_2e2m_wos=new double[nReviewPoints];

	double *signi_mean = new double[nReviewPoints];
	double *signi_median = new double[nReviewPoints];
	double *signi_m2s= new double[nReviewPoints];
	double *signi_m1s = new double[nReviewPoints];
	double *signi_p1s = new double[nReviewPoints];
	double *signi_p2s = new double[nReviewPoints];

	double *signi_mean_wos = new double[nReviewPoints];
	double *signi_median_wos = new double[nReviewPoints];
	double *signi_m2s_wos = new double[nReviewPoints];
	double *signi_m1s_wos = new double[nReviewPoints];
	double *signi_p1s_wos = new double[nReviewPoints];
	double *signi_p2s_wos = new double[nReviewPoints];

	double *shape_signi_mean_sys = new double[nReviewPoints];
	double *shape_signi_median_sys = new double[nReviewPoints];
	double *shape_cls_mean_sys = new double[nReviewPoints];
	double *shape_cls_r95_mean_sys = new double[nReviewPoints];

	double *shape_signi_mean_nosys = new double[nReviewPoints];
	double *shape_signi_median_nosys = new double[nReviewPoints];
	double *shape_cls_mean_nosys = new double[nReviewPoints];
	double *shape_cls_r95_mean_nosys = new double[nReviewPoints];

	double shape_bys_r95_sys[5][nReviewPoints];
	double shape_cls_r95_sys[5][nReviewPoints];
	double shape_cls_r95_nosys[5][nReviewPoints];
	double shape_bys_r95_nosys[5][nReviewPoints];
	double shape_signi_sys[6][nReviewPoints];
	double shape_signi_nosys[6][nReviewPoints];

	TString shellcommand;
	for(int i=0; i<nReviewPoints; i++){
		TString s; 
		for(int nsigma=0; nsigma<5; nsigma++){
			sprintf(ctmp, "grep CombinedR %s%1d | grep CLs_R95 | awk '{print $%d}'", s_log_nosys.c_str(), i, nsigma+11);
			s= gSystem->GetFromPipe(ctmp); 
			tmp_cls_r95_nosys[i][nsigma]=s.Atof();

			sprintf(ctmp, "grep CombinedR %s%1d | grep CLs_R95 | awk '{print $%d}'", s_log_fullysys.c_str(), i, nsigma+11);
			s= gSystem->GetFromPipe(ctmp); 
			tmp_cls_r95[i][nsigma]=s.Atof();

			sprintf(ctmp, "grep CombinedR %s%1d | grep Bys_R95 | awk '{print $%d}'", s_log_nosys.c_str(), i, nsigma+11);
			s= gSystem->GetFromPipe(ctmp); 
			tmp_bys_r95_nosys[i][nsigma]=s.Atof();

			sprintf(ctmp, "grep CombinedR %s%1d | grep Bys_R95 | awk '{print $%d}'", s_log_fullysys.c_str(), i, nsigma+11);
			s= gSystem->GetFromPipe(ctmp); 
			tmp_bys_r95[i][nsigma]=s.Atof();

			sprintf(ctmp, "grep CombinedR %s%1d | grep CLs_m2lnQ_b | awk '{print $%d}'", s_log_fullsys_correctedM2lnQ.c_str(), i, nsigma+4);
			s= gSystem->GetFromPipe(ctmp); 
			tmp_m2lnQ_b[i][nsigma]=s.Atof();

			sprintf(ctmp, "grep CombinedR %s%1d | grep CLs_m2lnQ_sb | awk '{print $%d}'", s_log_fullsys_correctedM2lnQ.c_str(), i, nsigma+4);
			s= gSystem->GetFromPipe(ctmp); 
			tmp_m2lnQ_sb[i][nsigma]=s.Atof();

			sprintf(ctmp, "grep CombinedR %s%1d | grep CLs_CLs | awk '{print $%d}'", s_log_fullysys.c_str(), i, nsigma+4);
			s= gSystem->GetFromPipe(ctmp); 
			tmp_cls_cls[i][nsigma]=s.Atof();
		}
		sprintf(ctmp, "grep CombinedR %s%1d | grep Bys_R95 | awk '{print $%d}'", s_log_fullysys.c_str(), i, 5+4);
		s= gSystem->GetFromPipe(ctmp); 
		tmp_bys_r95[i][5]=s.Atof();

		sprintf(ctmp, "grep CombinedR %s%1d | grep CLs_CLs | awk '{print $%d}'", s_log_4e.c_str(), i, 6);
		s= gSystem->GetFromPipe(ctmp); 
		cls_4e[i]=s.Atof();
		sprintf(ctmp, "grep CombinedR %s%1d | grep CLs_CLs | awk '{print $%d}'", s_log_4m.c_str(), i, 6);
		s= gSystem->GetFromPipe(ctmp); 
		cls_4m[i]=s.Atof();
		sprintf(ctmp, "grep CombinedR %s%1d | grep CLs_CLs | awk '{print $%d}'", s_log_2e2m.c_str(), i, 6);
		s= gSystem->GetFromPipe(ctmp); 
		cls_2e2m[i]=s.Atof();

		sprintf(ctmp, "grep CombinedR %s%1d | grep CLs_R95 | awk '{print $%d}'", s_log_4e.c_str(), i, 2+11);
		s= gSystem->GetFromPipe(ctmp); 
		clsr95_4e[i]=s.Atof();
		sprintf(ctmp, "grep CombinedR %s%1d | grep CLs_R95 | awk '{print $%d}'", s_log_4m.c_str(), i, 2+11);
		s= gSystem->GetFromPipe(ctmp); 
		clsr95_4m[i]=s.Atof();
		sprintf(ctmp, "grep CombinedR %s%1d | grep CLs_R95 | awk '{print $%d}'", s_log_2e2m.c_str(), i, 2+11);
		s= gSystem->GetFromPipe(ctmp); 
		clsr95_2e2m[i]=s.Atof();

		sprintf(ctmp, "grep CombinedR %s%1d | grep Bys_R95 | awk '{print $%d}'", s_log_4e_nosys.c_str(), i, 2+11);
		s= gSystem->GetFromPipe(ctmp); 
		bysr95_4e_wos[i]=s.Atof();
		sprintf(ctmp, "grep CombinedR %s%1d | grep Bys_R95 | awk '{print $%d}'", s_log_4m_nosys.c_str(), i, 2+11);
		s= gSystem->GetFromPipe(ctmp); 
		bysr95_4m_wos[i]=s.Atof();
		sprintf(ctmp, "grep CombinedR %s%1d | grep Bys_R95 | awk '{print $%d}'", s_log_2e2m_nosys.c_str(), i, 2+11);
		s= gSystem->GetFromPipe(ctmp); 
		bysr95_2e2m_wos[i]=s.Atof();

		sprintf(ctmp, "grep signi1fb_%d %s | awk '{print $4}'",i, s_log_signi.Data());
		s= gSystem->GetFromPipe(ctmp); 
		signi_m2s[i]=s.Atof();
		sprintf(ctmp, "grep signi1fb_%d %s | awk '{print $5}'",i, s_log_signi.Data());
		s= gSystem->GetFromPipe(ctmp); 
		signi_m1s[i]=s.Atof();
		sprintf(ctmp, "grep signi1fb_%d %s | awk '{print $7}'",i, s_log_signi.Data());
		s= gSystem->GetFromPipe(ctmp); 
		signi_p1s[i]=s.Atof();
		sprintf(ctmp, "grep signi1fb_%d %s | awk '{print $8}'",i, s_log_signi.Data());
		s= gSystem->GetFromPipe(ctmp); 
		signi_p2s[i]=s.Atof();
		sprintf(ctmp, "grep signi1fb_%d %s | awk '{print $6}'",i, s_log_signi.Data());
		s= gSystem->GetFromPipe(ctmp); 
		signi_median[i]=s.Atof();
		sprintf(ctmp, "grep signi1fb_%d %s | awk '{print $10}'",i, s_log_signi.Data());
		s= gSystem->GetFromPipe(ctmp); 
		signi_mean[i]=s.Atof();

		sprintf(ctmp, "grep signi_nosys_%d %s | awk '{print $4}'",i, s_log_signi_nosys.Data());
		s= gSystem->GetFromPipe(ctmp); 
		signi_m2s_wos[i]=s.Atof();
		sprintf(ctmp, "grep signi_nosys_%d %s | awk '{print $5}'",i, s_log_signi_nosys.Data());
		s= gSystem->GetFromPipe(ctmp); 
		signi_m1s_wos[i]=s.Atof();
		sprintf(ctmp, "grep signi_nosys_%d %s | awk '{print $7}'",i, s_log_signi_nosys.Data());
		s= gSystem->GetFromPipe(ctmp); 
		signi_p1s_wos[i]=s.Atof();
		sprintf(ctmp, "grep signi_nosys_%d %s | awk '{print $8}'",i, s_log_signi_nosys.Data());
		s= gSystem->GetFromPipe(ctmp); 
		signi_p2s_wos[i]=s.Atof();
		sprintf(ctmp, "grep signi_nosys_%d %s | awk '{print $6}'",i, s_log_signi_nosys.Data());
		s= gSystem->GetFromPipe(ctmp); 
		signi_median_wos[i]=s.Atof();
		sprintf(ctmp, "grep signi_nosys_%d %s | awk '{print $10}'",i, s_log_signi_nosys.Data());
		s= gSystem->GetFromPipe(ctmp); 
		signi_mean_wos[i]=s.Atof();

		shellcommand.Form("grep ToS_significance %s | grep -E \"mH=%3.0f\" | awk '{print $6}' ", s_log_shape_sys.Data(), xReviewPoints[i]  );
		s=gSystem->GetFromPipe(shellcommand);
		shape_signi_median_sys[i]=s.Atof();
		shellcommand.Form("grep ToS_significance %s | grep -E \"mH=%3.0f\" | awk '{print $9}' ", s_log_shape_sys.Data(), xReviewPoints[i]  );
		s=gSystem->GetFromPipe(shellcommand);
		shape_signi_mean_sys[i]=s.Atof();
		shellcommand.Form("grep CLs_CLs %s | grep -E \"mH=%3.0f\" | awk '{print $6}' ", s_log_shape_sys.Data(), xReviewPoints[i]  );
		s=gSystem->GetFromPipe(shellcommand);
		shape_cls_mean_sys[i]=s.Atof();
		shellcommand.Form("grep CLs_R95 %s | grep -E \"mH=%3.0f\" | awk '{print $13}' ", s_log_shape_sys.Data(), xReviewPoints[i]  );
		s=gSystem->GetFromPipe(shellcommand);
		shape_cls_r95_mean_sys[i]=s.Atof();

		shellcommand.Form("grep ToS_significance %s | grep -E \"mH=%3.0f\" | awk '{print $6}' ", s_log_shape_nosys.Data(), xReviewPoints[i]  );
		s=gSystem->GetFromPipe(shellcommand);
		shape_signi_median_nosys[i]=s.Atof();
		shellcommand.Form("grep ToS_significance %s | grep -E \"mH=%3.0f\" | awk '{print $9}' ", s_log_shape_nosys.Data(), xReviewPoints[i]  );
		s=gSystem->GetFromPipe(shellcommand);
		shape_signi_mean_nosys[i]=s.Atof();
		shellcommand.Form("grep CLs_CLs %s | grep -E \"mH=%3.0f\" | awk '{print $6}' ", s_log_shape_nosys.Data(), xReviewPoints[i]  );
		s=gSystem->GetFromPipe(shellcommand);
		shape_cls_mean_nosys[i]=s.Atof();
		shellcommand.Form("grep CLs_R95 %s | grep -E \"mH=%3.0f\" | awk '{print $13}' ", s_log_shape_nosys.Data(), xReviewPoints[i]  );
		s=gSystem->GetFromPipe(shellcommand);
		shape_cls_r95_mean_nosys[i]=s.Atof();

		for(int sigma=0; sigma<6; sigma++){
			shellcommand.Form("grep ToS_significance %s | grep -E \"mH=%3.0f\" | awk '{print $%d}' ", s_log_shape_sys.Data(), xReviewPoints[i] , sigma+4 );
			s=gSystem->GetFromPipe(shellcommand);
			shape_signi_sys[sigma][i]=s.Atof();
			shellcommand.Form("grep ToS_significance %s | grep -E \"mH=%3.0f\" | awk '{print $%d}' ", s_log_shape_nosys.Data(), xReviewPoints[i] , sigma+4 );
			s=gSystem->GetFromPipe(shellcommand);
			shape_signi_nosys[sigma][i]=s.Atof();


			if(sigma>=5) continue;
			shellcommand.Form("grep CLs_R95 %s | grep -E \"mH=%3.0f\" | awk '{print $%d}' ", s_log_shape_sys.Data(), xReviewPoints[i], sigma+11  );
			s=gSystem->GetFromPipe(shellcommand);
			shape_cls_r95_sys[sigma][i]=s.Atof();
			shellcommand.Form("grep Bys_R95 %s | grep -E \"mH=%3.0f\" | awk '{print $%d}' ", s_log_shape_sys.Data(), xReviewPoints[i], sigma+11  );
			s=gSystem->GetFromPipe(shellcommand);
			shape_bys_r95_sys[sigma][i]=s.Atof();


			shellcommand.Form("grep CLs_R95 %s | grep -E \"mH=%3.0f\" | awk '{print $%d}' ", s_log_shape_nosys.Data(), xReviewPoints[i], sigma+11  ); // 11 for 200pb, 12 for 1fb
			s=gSystem->GetFromPipe(shellcommand);
			shape_cls_r95_nosys[sigma][i]=s.Atof();
			shellcommand.Form("grep Bys_R95 %s | grep -E \"mH=%3.0f\" | awk '{print $%d}' ", s_log_shape_nosys.Data(), xReviewPoints[i], sigma+11  );
			s=gSystem->GetFromPipe(shellcommand);
			shape_bys_r95_nosys[sigma][i]=s.Atof();
		}

	}

	cout<<"nReviewPoints"<<nReviewPoints<<endl;


	double *cmb_rBys_expB= new double[nReviewPoints];
	double *cmb_rBys_mean = new double[nReviewPoints];
	double *cmb_rFrq_mean = new double[nReviewPoints];
	double *cmb_rBys_m1s = new double[nReviewPoints];
	double *cmb_rFrq_m1s = new double[nReviewPoints];
	double *cmb_rBys_m2s = new double[nReviewPoints];
	double *cmb_rFrq_m2s = new double[nReviewPoints];
	double *cmb_rBys_p1s = new double[nReviewPoints];
	double *cmb_rFrq_p1s = new double[nReviewPoints];
	double *cmb_rBys_p2s = new double[nReviewPoints];
	double *cmb_rFrq_p2s = new double[nReviewPoints];

	double *cmb_m2lnQ_b_sigma[5], *cmb_m2lnQ_sb_sigma[5];//index from 0 to 4:  -2s, -1s, mean, 1s, 2s
	double *cmb_cls_sigma[5];//index from 0 to 4:  -2s, -1s, mean, 1s, 2s   // you should notice that it's mean not median
	double *cmb_cls_r95_nosys_sigma[5];
	double *cmb_bys_r95_nosys_sigma[5];
	for(int i=0; i<5; i++){
		cmb_m2lnQ_b_sigma[i]=new double[nReviewPoints];
		cmb_m2lnQ_sb_sigma[i]=new double[nReviewPoints];
		cmb_cls_sigma[i]=new double[nReviewPoints];
		cmb_cls_r95_nosys_sigma[i]=new double[nReviewPoints];
		cmb_bys_r95_nosys_sigma[i]=new double[nReviewPoints];
	}

	for(int i=0; i<nReviewPoints; i++){
		cmb_rBys_expB[i]=tmp_bys_r95[i][5];
		cmb_rBys_mean[i]=tmp_bys_r95[i][2];
		cmb_rBys_m1s[i]=tmp_bys_r95[i][1];
		cmb_rBys_m2s[i]=tmp_bys_r95[i][0];
		cmb_rBys_p1s[i]=tmp_bys_r95[i][3];
		cmb_rBys_p2s[i]=tmp_bys_r95[i][4];
		cmb_rFrq_mean[i]=tmp_cls_r95[i][2];
		cmb_rFrq_m1s[i]=tmp_cls_r95[i][1];
		cmb_rFrq_m2s[i]=tmp_cls_r95[i][0];
		cmb_rFrq_p1s[i]=tmp_cls_r95[i][3];
		cmb_rFrq_p2s[i]=tmp_cls_r95[i][4];
		for(int sigma=0; sigma<5; sigma++){
			cmb_m2lnQ_b_sigma[sigma][i]=tmp_m2lnQ_b[i][sigma];
			cmb_m2lnQ_sb_sigma[sigma][i]=tmp_m2lnQ_sb[i][sigma];
			cmb_cls_sigma[sigma][i]=tmp_cls_cls[i][sigma];
			cmb_cls_r95_nosys_sigma[sigma][i]=tmp_cls_r95_nosys[i][sigma];
			cmb_bys_r95_nosys_sigma[sigma][i]=tmp_bys_r95_nosys[i][sigma];
		}
	}


	TPaveText *pt;

	for(int np=0; np<nReviewPoints; np++) {
		if(bProcessMultiChannels && debug ){
			cout<<" Combined Channels: "<<endl;
			double rm1s=0, rm2s=0, rp1s=0, rp2s=0, r0s=0, rmean=0, rExpB=0; 
			rmean=cmb_rFrq_mean[np]; 
			rm1s=cmb_rFrq_m1s[np]; 
			rm2s=cmb_rFrq_m2s[np]; 
			rp1s=cmb_rFrq_p1s[np]; 
			rp2s=cmb_rFrq_p2s[np]; 

			printf("CombinedResultsCLs_R95 -- mH=%3.0fGeV %10.5f %10.5f %10.5f %10.5f %10.5f\n",
					xReviewPoints[np],rm2s,rm1s,
					rmean, rp1s,rp2s);

			rmean=cmb_cls_sigma[2][np];
			rm2s=cmb_cls_sigma[0][np];
			rm1s=cmb_cls_sigma[1][np];
			rp1s=cmb_cls_sigma[3][np];
			rp2s=cmb_cls_sigma[4][np];

			printf("CombinedResultsCLs_CLs -- mH=%3.0fGeV %10.5f %10.5f %10.5f %10.5f %10.5f\n",
					xReviewPoints[np], rm2s, rm1s, rmean, rp1s, rp2s);

			rmean=cmb_m2lnQ_b_sigma[2][np];
			rm2s=cmb_m2lnQ_b_sigma[0][np];
			rm1s=cmb_m2lnQ_b_sigma[1][np];
			rp1s=cmb_m2lnQ_b_sigma[3][np];
			rp2s=cmb_m2lnQ_b_sigma[4][np];
			if(debug)printf("CombinedResultsCLs_m2lnQ_b -- mH=%3.0fGeV %10.5f %10.5f %10.5f %10.5f %10.5f\n",
					xReviewPoints[np], rm2s, rm1s, rmean, rp1s, rp2s);

			rmean=cmb_m2lnQ_sb_sigma[2][np];
			rm2s=cmb_m2lnQ_sb_sigma[0][np];
			rm1s=cmb_m2lnQ_sb_sigma[1][np];
			rp1s=cmb_m2lnQ_sb_sigma[3][np];
			rp2s=cmb_m2lnQ_sb_sigma[4][np];
			printf("CombinedResultsCLs_m2lnQ_sb -- mH=%3.0fGeV %10.5f %10.5f %10.5f %10.5f %10.5f\n",
					xReviewPoints[np], rm2s, rm1s, rmean, rp1s, rp2s);

			rm1s=cmb_rBys_m1s[np];
			rm2s=cmb_rBys_m2s[np];
			rp1s=cmb_rBys_p1s[np];
			rp2s=cmb_rBys_p2s[np];
			rmean=cmb_rBys_mean[np];
			rExpB =cmb_rBys_expB[np];
			printf("CombinedResultsBys_R95 -- mH=%3.0fGeV %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f\n",
					xReviewPoints[np],rm2s,rm1s,
					rmean, rp1s,rp2s, rExpB);

		}//end bProcessMultiChannels

	}//end loop points

	//--------------------------show the comparison of exclusion limits results from Bayesian and CLs as function of Higgs mass
	if(bProcessMultiChannels) {


		//----for combined plots........
		if(s_text_topleft.Contains("200 pb") && projectingRLimitLogY ) 
			pt = SetTPaveText(0.2, 0.55, 0.55, 0.60);
		else 
			pt = SetTPaveText(0.5, 0.7, 0.8, 0.9);
		pt->AddText(schannels.c_str());
		sprintf(ctmp,"%s_count_combined_cls",projectingRLimitSavePath.Data());
		ssave=ctmp;
		PlotWithBelts combinedlimits2(
				cmb_rFrq_m1s, cmb_rFrq_p1s, cmb_rFrq_m2s, cmb_rFrq_p2s,
				//cmb_rFrq_mean, cmb_rBys_mean, nReviewPoints,  
				cmb_rFrq_mean, cmb_cls_r95_nosys_sigma[2], nReviewPoints,  
				xReviewPoints, ssave, pt,
				projectingXmin, projectingXmax, projectingRLimitYmin, projectingRLimitYmax, projectingRLimitLogY, projectingRLimitXYtitles);
		combinedlimits2.plot();
		combinedlimits2.getObsGraph()->SetMarkerStyle(1);//remove marker
		//combinedlimits2.drawLegend("CLs bkgd-only: mean","CLs bkgd-only: 1 #sigma band", "CLs bkgd-only: 2 #sigma band", "Bayesian bkgd-only: mean (no sys)");
		combinedlimits2.drawLegend("CLs 95% CL exclusion: mean","95% CL exclusion: 68% band", "95% CL exclusion: 95% band", "95% CL exclusion: mean (no sys)");
		MoveLegend(combinedlimits2.getLegend(),0.5,0.6);
		if(s_text_topleft.Contains("200 pb") && projectingRLimitLogY ) MoveLegend(combinedlimits2.getLegend(),0.2,0.35);
		combinedlimits2.getMeanGraph()->SetLineStyle(1);
		combinedlimits2.getMeanGraph()->SetLineColor(kBlue);
		combinedlimits2.getMeanGraph()->SetLineWidth(2);
		combinedlimits2.getObsGraph()->SetLineStyle(2);
		combinedlimits2.getObsGraph()->SetLineWidth(1);
		combinedlimits2.getLine()->SetLineColor(kRed);
		DrawText(0.145, 0.9, 0.8, 1.0, s_text_topleft);
		DrawText(0.82, 0.9, 1.0, 1.0, s_text_topright, 0.04);
		combinedlimits2.save();	

		//----for combined plots........
		if(s_text_topleft.Contains("200 pb") && projectingRLimitLogY ) 
			pt = SetTPaveText(0.2, 0.55, 0.55, 0.60);
		else 
			pt = SetTPaveText(0.5, 0.7, 0.8, 0.9);
		pt->AddText(schannels.c_str());
		sprintf(ctmp,"%s_count_combined_bayesian",projectingRLimitSavePath.Data());
		ssave=ctmp;
		PlotWithBelts combinedlimits3(
				cmb_rBys_m1s, cmb_rBys_p1s, cmb_rBys_m2s, cmb_rBys_p2s,
				//cmb_rBys_mean, cmb_rBys_mean, nReviewPoints,  
				cmb_rBys_mean, cmb_bys_r95_nosys_sigma[2], nReviewPoints,  
				xReviewPoints, ssave, pt,
				projectingXmin, projectingXmax, projectingRLimitYmin, projectingRLimitYmax, projectingRLimitLogY, projectingRLimitXYtitles);
		combinedlimits3.plot();
		combinedlimits3.getObsGraph()->SetMarkerStyle(1);//remove marker
		//combinedlimits3.drawLegend("CLs bkgd-only: mean","CLs bkgd-only: 1 #sigma band", "CLs bkgd-only: 2 #sigma band", "Bayesian bkgd-only: mean (no sys)");
		combinedlimits3.drawLegend("Bayesian 95% CL exclusion: mean","95% CL exclusion: 68% band", "95% CL exclusion: 95% band", "95% CL exclusion: mean (no sys)");
		MoveLegend(combinedlimits3.getLegend(),0.5,0.6);
		if(s_text_topleft.Contains("200 pb") && projectingRLimitLogY ) MoveLegend(combinedlimits3.getLegend(),0.2,0.35);
		combinedlimits3.getMeanGraph()->SetLineStyle(1);
		combinedlimits3.getMeanGraph()->SetLineColor(kBlue);
		combinedlimits3.getMeanGraph()->SetLineWidth(2);
		combinedlimits3.getObsGraph()->SetLineStyle(2);
		combinedlimits3.getObsGraph()->SetLineWidth(1);
		combinedlimits3.getLine()->SetLineColor(kRed);
		DrawText(0.145, 0.9, 0.8, 1.0, s_text_topleft);
		DrawText(0.82, 0.9, 1.0, 1.0, s_text_topright, 0.04);
		combinedlimits3.save();	

		pt = SetTPaveText(0.2, 0.35, 0.5, 0.45);
		pt->AddText("CLs (with sys)");
		pt->AddText(schannels.c_str());
		sprintf(ctmp,"%s_combined_count", projectingCLsSavePath.Data());
		ssave=ctmp;
		PlotWithBelts clswithbets(cmb_cls_sigma[1], cmb_cls_sigma[3], cmb_cls_sigma[0], cmb_cls_sigma[4],
				cmb_cls_sigma[2], nReviewPoints, xReviewPoints, ssave, pt,
				projectingXmin, projectingXmax, projectingCLsYmin, projectingCLsYmax, projectingCLsLogY, projectingCLsXYtitles);
		clswithbets.plot();
		clswithbets.drawLegend("bkg-only: mean","bkg-only: 68% band", "bkg-only: 95% band");
		MoveLegend(clswithbets.getLegend(),0.2,0.2);
		DrawText(0.145, 0.9, 0.8, 1.0, s_text_topleft);
		DrawText(0.82, 0.9, 1.0, 1.0, s_text_topright, 0.04);
		clswithbets.save();

		pt = SetTPaveText(0.2, 0.35, 0.5, 0.45);
		pt->AddText("CLs (with sys)");
		pt->AddText(schannels.c_str());
		sprintf(ctmp,"%s_eachchannel+combined_count",projectingCLsSavePath.Data());
		ssave=ctmp;
		DrawMultiGraph dmgclsseparate(
				xReviewPoints, cmb_cls_sigma[2], nReviewPoints, "Combination",
				xReviewPoints, cls_4e, nReviewPoints, "4e",
				projectingCLsXYtitles, ssave, pt, 
				projectingXmin, projectingXmax, 
				0.03, 1, 1);
		dmgclsseparate.add(xReviewPoints, cls_4m, nReviewPoints, "4m");
		dmgclsseparate.add(xReviewPoints, cls_2e2m, nReviewPoints, "2e2m");
		dmgclsseparate.draw();
		DrawText(0.145, 0.9, 0.8, 1.0, s_text_topleft);
		DrawText(0.82, 0.9, 1.0, 1.0, s_text_topright, 0.04);
		MoveLegend(dmgclsseparate.getLegend(), 0.2, 0.2);
		Save(dmgclsseparate.getCanvas(), ssave);
/*
		clswithbets.getCanvas()->cd();
		dmgclsseparate.getGraph(1)->Draw("sameLP");
		dmgclsseparate.getGraph(2)->Draw("sameLP");
		dmgclsseparate.getGraph(3)->Draw("sameLP");
		dmgclsseparate.getLegend()->Draw();
		Save(clswithbets.getCanvas(), ssave);
*/
		pt = SetTPaveText(0.67, 0.7, 0.8, 0.9);
		pt->AddText("-2lnQ (with sys)");
		pt->AddText(schannels.c_str());
		sprintf(ctmp,"%s_b_combined_count", projectingM2lnQSavePath.Data());
		ssave=ctmp;
		PlotWithBelts m2lnQwithbelts_b(
				cmb_m2lnQ_b_sigma[1], cmb_m2lnQ_b_sigma[3], cmb_m2lnQ_b_sigma[0], cmb_m2lnQ_b_sigma[4],
				cmb_m2lnQ_b_sigma[2], cmb_m2lnQ_sb_sigma[2],nReviewPoints,
				xReviewPoints, ssave, pt,
				projectingXmin, projectingXmax, projectingM2lnQYmin, projectingM2lnQYmax, projectingM2lnQLogY,
				projectingM2lnQXYtitles);
		m2lnQwithbelts_b.plot();
		m2lnQwithbelts_b.getObsGraph()->SetMarkerStyle(1);//no marker
		m2lnQwithbelts_b.drawLegend("bkgd-only: mean","bkgd-only: 68% band", "bkgd-only: 95% band", "bkgd+signal: mean");
		m2lnQwithbelts_b.getLine()->SetY1(0.);
		m2lnQwithbelts_b.getLine()->SetY2(0.);
		MoveLegend(m2lnQwithbelts_b.getLegend(),0.5,0.3);
		m2lnQwithbelts_b.getLine()->SetLineColor(kRed);
		DrawText(0.145, 0.9, 0.8, 1.0, s_text_topleft);
		DrawText(0.82, 0.9, 1.0, 1.0, s_text_topright, 0.04);
		m2lnQwithbelts_b.save();	

		sprintf(ctmp,"%s_sb_combined_count", projectingM2lnQSavePath.Data());
		ssave=ctmp;
		PlotWithBelts m2lnQwithbelts_sb(
				cmb_m2lnQ_sb_sigma[1], cmb_m2lnQ_sb_sigma[3], cmb_m2lnQ_sb_sigma[0], cmb_m2lnQ_sb_sigma[4],
				cmb_m2lnQ_b_sigma[2], cmb_m2lnQ_sb_sigma[2],nReviewPoints,
				xReviewPoints, ssave, pt,
				projectingXmin, projectingXmax, projectingM2lnQYmin, projectingM2lnQYmax, projectingM2lnQLogY,
				projectingM2lnQXYtitles);	
		m2lnQwithbelts_sb.plot();
		m2lnQwithbelts_sb.getObsGraph()->SetMarkerStyle(1);//no marker
		m2lnQwithbelts_sb.drawLegend("bkgd-only: mean","bkgd+signal: 68% band", "bkgd+signal: 95% band", "bkgd+signal: mean");
		m2lnQwithbelts_sb.getLine()->SetY1(0.);
		m2lnQwithbelts_sb.getLine()->SetY2(0.);
		MoveLegend(m2lnQwithbelts_sb.getLegend(),0.2,0.6);
		m2lnQwithbelts_sb.getLine()->SetLineColor(kRed);
		DrawText(0.145, 0.9, 0.8, 1.0, s_text_topleft);
		DrawText(0.82, 0.9, 1.0, 1.0, s_text_topright, 0.04);
		m2lnQwithbelts_sb.save();	

		pt = SetTPaveText(0.7, 0.8, 0.8, 0.9);
		pt->AddText(schannels.c_str());
		sprintf(ctmp,"%s_cls_eachchannel+combined_count",projectingRLimitSavePath.Data());
		ssave=ctmp;
		DrawMultiGraph dmg(
				xReviewPoints, cmb_rFrq_mean, nReviewPoints, "Combination",
				xReviewPoints, clsr95_4e, nReviewPoints, "4e",
				projectingRLimitXYtitles, ssave, pt, 
				projectingXmin, projectingXmax, 
				0.5, 100, 1);
		dmg.add(xReviewPoints, clsr95_4m, nReviewPoints, "4m");
		dmg.add(xReviewPoints, clsr95_2e2m, nReviewPoints, "2e2m");
		dmg.draw();
		DrawText(0.145, 0.9, 0.8, 1.0, s_text_topleft);
		DrawText(0.82, 0.9, 1.0, 1.0, s_text_topright, 0.04);
		MoveLegend(dmg.getLegend(), 0.7, 0.65);
		Save(dmg.getCanvas(), ssave);

		pt = SetTPaveText(0.7, 0.8, 0.8, 0.9);
		pt->AddText(schannels.c_str());
		sprintf(ctmp,"%s_bayesian_eachchannel+combined_count",projectingRLimitSavePath.Data());
		ssave=ctmp;
		DrawMultiGraph dmgBayesianEachComLimit(
				xReviewPoints, cmb_bys_r95_nosys_sigma[2], nReviewPoints, "Combination (no sys)",
				xReviewPoints, bysr95_4e_wos, nReviewPoints, "4e (no sys)",
				projectingRLimitXYtitles, ssave, pt, 
				projectingXmin, projectingXmax, 
				0.5, 100, 1);
		dmgBayesianEachComLimit.add(xReviewPoints, bysr95_4m_wos, nReviewPoints, "4m (no sys)");
		dmgBayesianEachComLimit.add(xReviewPoints, bysr95_2e2m_wos, nReviewPoints, "2e2m (no sys)");
		dmgBayesianEachComLimit.draw();
		DrawText(0.145, 0.9, 0.8, 1.0, s_text_topleft);
		DrawText(0.82, 0.9, 1.0, 1.0, s_text_topright, 0.04);
		MoveLegend(dmgBayesianEachComLimit.getLegend(), 0.7, 0.65);
		Save(dmgBayesianEachComLimit.getCanvas(), ssave);

		if(s_text_topleft.Contains("1 fb") )  {
			pt = SetTPaveText(0.7, 0.8, 0.8, 0.9);
			pt->AddText(schannels.c_str());
			sprintf(ctmp,"%s_mean_count_combined_SYS.vs.NOSYS",projectingSignificanceSavePath.Data());
			ssave=ctmp;
			DrawMultiGraph dmgsigni(
					xReviewPoints, signi_mean, nReviewPoints, "mean (w/ sys)",
					xReviewPoints, signi_mean_wos, nReviewPoints, "mean (no sys)",
					projectingSignificanceXYtitles, ssave, pt, 
					projectingXmin, projectingXmax, 
					projectingSignificanceYmin, projectingSignificanceYmax, projectingSignificanceLogY);
			dmgsigni.draw();
			DrawText(0.145, 0.9, 0.8, 1.0, s_text_topleft);
			DrawText(0.82, 0.9, 1.0, 1.0, s_text_topright, 0.04);
			MoveLegend(dmgsigni.getLegend(), 0.7, 0.65);
			dmgsigni.getLine()->Delete();
			Save(dmgsigni.getCanvas(), ssave);

			pt = SetTPaveText(0.4, 0.85, 0.8, 0.9);
			pt->AddText(schannels.c_str());
			sprintf(ctmp,"%s_count_combined_bands",projectingSignificanceSavePath.Data());
			ssave=ctmp;
			PlotWithBelts beltSigni(
					signi_m1s, signi_p1s, signi_m2s, signi_p2s,
					signi_mean, signi_mean_wos, nReviewPoints,  
					xReviewPoints, ssave, pt,
					projectingXmin, projectingXmax, projectingSignificanceYmin, projectingSignificanceYmax, projectingSignificanceLogY, projectingSignificanceXYtitles);
			beltSigni.plot();
			beltSigni.getObsGraph()->SetMarkerStyle(1);//remove marker
			beltSigni.drawLegend("mean","68% band", "95% band", "mean(no sys)");
			MoveLegend(beltSigni.getLegend(),0.4,0.7);
			//if(s_text_topleft.Contains("200 pb") && projectingSignificanceLogY ) MoveLegend(beltSigni.getLegend(),0.2,0.35);
			beltSigni.getMeanGraph()->SetLineStyle(1);
			beltSigni.getMeanGraph()->SetLineColor(kBlue);
			beltSigni.getMeanGraph()->SetLineWidth(2);
			beltSigni.getObsGraph()->SetLineStyle(2);
			beltSigni.getObsGraph()->SetLineWidth(1);
			beltSigni.getLine()->SetLineColor(kRed);
			DrawText(0.145, 0.9, 0.8, 1.0, s_text_topleft);
			DrawText(0.82, 0.9, 1.0, 1.0, s_text_topright, 0.04);
			beltSigni.getLine()->Delete();
			beltSigni.save();	

			pt = SetTPaveText(0.7, 0.8, 0.8, 0.9);
			pt->AddText(schannels.c_str());
			sprintf(ctmp,"%s_mean_Count.vs.Shape_withSys",projectingSignificanceSavePath.Data());
			ssave=ctmp;
			DrawMultiGraph dmgsigCSmean(
					xReviewPoints, shape_signi_mean_sys, nReviewPoints, "mean (shape w/ sys)",
					xReviewPoints, signi_mean, nReviewPoints, "mean (count w/ sys)",
					projectingSignificanceXYtitles, ssave, pt, 
					projectingXmin, projectingXmax, 
					projectingSignificanceYmin, projectingSignificanceYmax, projectingSignificanceLogY);
			dmgsigCSmean.draw();
			DrawText(0.145, 0.9, 0.8, 1.0, s_text_topleft);
			DrawText(0.82, 0.9, 1.0, 1.0, s_text_topright, 0.04);
			MoveLegend(dmgsigCSmean.getLegend(), 0.5, 0.65);
			dmgsigCSmean.getLine()->Delete();
			Save(dmgsigCSmean.getCanvas(), ssave);
		}
		pt = SetTPaveText(0.5, 0.4, 0.8, 0.6);
		pt->AddText(schannels.c_str());
		sprintf(ctmp,"%s_mean_Count.vs.Shape_withSys", projectingCLsSavePath.Data());
		ssave=ctmp;
		DrawMultiGraph dmgclsCSmean(
				xReviewPoints, shape_cls_mean_sys, nReviewPoints, "mean (shape w/ sys)",
				xReviewPoints, cmb_cls_sigma[2], nReviewPoints, "mean (count w/ sys)",
				projectingCLsXYtitles,ssave, pt,
				projectingXmin, projectingXmax, projectingCLsYmin, projectingCLsYmax, projectingCLsLogY);
		dmgclsCSmean.draw();
		MoveLegend(dmgclsCSmean.getLegend(),0.5,0.3);
		DrawText(0.145, 0.9, 0.8, 1.0, s_text_topleft);
		DrawText(0.82, 0.9, 1.0, 1.0, s_text_topright, 0.04);
		dmgclsCSmean.getLine()->SetY1(0.05);
		dmgclsCSmean.getLine()->SetY2(0.05);
		Save(dmgclsCSmean.getCanvas(), ssave);

		pt = SetTPaveText(0.7, 0.8, 0.8, 0.9);
		pt->AddText(schannels.c_str());
		sprintf(ctmp,"%s_clsr95mean_Count.vs.Shape_withSys",projectingRLimitSavePath.Data());
		ssave=ctmp;
		DrawMultiGraph dmgclsr95CSmean(
				//xReviewPoints, shape_cls_r95_mean_sys, nReviewPoints, "mean (shape w/ sys)",
				//xReviewPoints, cmb_rFrq_mean, nReviewPoints, "mean (count w/ sys)",
				xReviewPoints, shape_cls_r95_mean_sys, nReviewPoints, "mean (30% sys on Higgs XS)",
				xReviewPoints, cmb_rFrq_mean, nReviewPoints, "mean (10% sys on HiggsXS)",
				projectingRLimitXYtitles, ssave, pt, 
				projectingXmin, projectingXmax, 
				0.5, 100, 1);
		dmgclsr95CSmean.draw();
		DrawText(0.145, 0.9, 0.8, 1.0, s_text_topleft);
		DrawText(0.82, 0.9, 1.0, 1.0, s_text_topright, 0.04);
		MoveLegend(dmgclsr95CSmean.getLegend(), 0.5, 0.65);
		Save(dmgclsr95CSmean.getCanvas(), ssave);

		pt = SetTPaveText(0.7, 0.8, 0.8, 0.9);
		pt->AddText(schannels.c_str());
		sprintf(ctmp,"%s_count_combined_Bayesian.vs.CLs_nosys",projectingRLimitSavePath.Data());
		ssave=ctmp;
		DrawMultiGraph dmgclsvsbys(
				xReviewPoints, cmb_bys_r95_nosys_sigma[2] , nReviewPoints, "Bayesian (no sys)",
				xReviewPoints, cmb_cls_r95_nosys_sigma[2] , nReviewPoints, "Modified Frequentist (no sys)",
				projectingRLimitXYtitles, ssave, pt, 
				projectingXmin, projectingXmax, 
				0.5, 100, 1);
		dmgclsvsbys.draw();
		DrawText(0.145, 0.9, 0.8, 1.0, s_text_topleft);
		DrawText(0.82, 0.9, 1.0, 1.0, s_text_topright, 0.04);
		MoveLegend(dmgclsvsbys.getLegend(), 0.5, 0.65);
		Save(dmgclsvsbys.getCanvas(), ssave);

		cout<<"\n\n";
		cout<<"Combined result: significance "<<endl;
		cout<<"        mass(GeV)      -2sigma     -1sigma     mean       1sigma       2sigma   median "<<endl;
		for(int np=0; np<nReviewPoints; np++){
			printf("w/ sys  %13.0f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f \n",
					xReviewPoints[np],
					signi_m2s[np],
					signi_m1s[np],
					signi_mean[np],
					signi_p1s[np],
					signi_p2s[np],
					signi_median[np]
			      );
			printf("no sys  %13.0f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f \n",
					xReviewPoints[np],
					signi_m2s_wos[np],
					signi_m1s_wos[np],
					signi_mean_wos[np],
					signi_p1s_wos[np],
					signi_p2s_wos[np],
					signi_median_wos[np]
			      );
		}

		// print out  in good way
		cout<<"\n\n";
		cout<<"Combined result: C.L. 95% limit "<<endl;
		cout<<"method    mass(GeV) -2sigma    -1sigma    mean       1sigma    2sigma "<<endl;
		for(int np=0; np<nReviewPoints; np++){
			printf("CLs(w/sys)   %3.0f %10.5f %10.5f %10.5f %10.5f %10.5f\n",
					xReviewPoints[np],
					cmb_rFrq_m2s[np],
					cmb_rFrq_m1s[np],
					cmb_rFrq_mean[np],
					cmb_rFrq_p1s[np],
					cmb_rFrq_p2s[np]);
			printf("CLs(nosys)   %3.0f %10.5f %10.5f %10.5f %10.5f %10.5f\n",
					xReviewPoints[np],
					cmb_cls_r95_nosys_sigma[0][np],
					cmb_cls_r95_nosys_sigma[1][np],
					cmb_cls_r95_nosys_sigma[2][np],
					cmb_cls_r95_nosys_sigma[3][np],
					cmb_cls_r95_nosys_sigma[4][np]
			      );
			printf("Bys(w/sys)   %3.0f %10.5f %10.5f %10.5f %10.5f %10.5f\n",
					xReviewPoints[np],
					cmb_rBys_m2s[np],
					cmb_rBys_m1s[np],
					cmb_rBys_mean[np],
					cmb_rBys_p1s[np],
					cmb_rBys_p2s[np]);
			printf("Bys(nosys)   %3.0f %10.5f %10.5f %10.5f %10.5f %10.5f\n",
					xReviewPoints[np],
					cmb_bys_r95_nosys_sigma[0][np],
					cmb_bys_r95_nosys_sigma[1][np],
					cmb_bys_r95_nosys_sigma[2][np],
					cmb_bys_r95_nosys_sigma[3][np],
					cmb_bys_r95_nosys_sigma[4][np]
			      );
		}
		cout<<"\n\n"<<endl;

		cout<<"\n\nSHAPE\n\n";
		cout<<"Combined result: significance "<<endl;
		cout<<"        mass(GeV)      -2sigma     -1sigma     mean       1sigma       2sigma   median "<<endl;
		for(int np=0; np<nReviewPoints; np++){
			printf("w/ sys  %13.0f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f \n",
					xReviewPoints[np],
					shape_signi_sys[0][np],
					shape_signi_sys[1][np],
					shape_signi_sys[5][np],
					shape_signi_sys[3][np],
					shape_signi_sys[4][np],
					shape_signi_sys[2][np]
			      );
			printf("no sys  %13.0f %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f \n",
					xReviewPoints[np],
					shape_signi_nosys[0][np],
					shape_signi_nosys[1][np],
					shape_signi_nosys[5][np],
					shape_signi_nosys[3][np],
					shape_signi_nosys[4][np],
					shape_signi_nosys[2][np]
			      );
		}

		// print out  in good way
		cout<<"\n\n";
		cout<<"Combined result: C.L. 95% limit "<<endl;
		cout<<"method    mass(GeV) -2sigma    -1sigma    mean       1sigma    2sigma "<<endl;
		for(int np=0; np<nReviewPoints; np++){
			printf("CLs(w/sys)   %3.0f %10.5f %10.5f %10.5f %10.5f %10.5f\n",
					xReviewPoints[np],
					shape_cls_r95_sys[0][np],
					shape_cls_r95_sys[1][np],
					shape_cls_r95_sys[2][np],
					shape_cls_r95_sys[3][np],
					shape_cls_r95_sys[4][np]
			      );
			printf("CLs(nosys)   %3.0f %10.5f %10.5f %10.5f %10.5f %10.5f\n",
					xReviewPoints[np],
					shape_cls_r95_nosys[0][np],
					shape_cls_r95_nosys[1][np],
					shape_cls_r95_nosys[2][np],
					shape_cls_r95_nosys[3][np],
					shape_cls_r95_nosys[4][np]
			      );
			printf("Bys(w/sys)   %3.0f %10.5f %10.5f %10.5f %10.5f %10.5f\n",
					xReviewPoints[np],
					shape_bys_r95_sys[0][np],
					shape_bys_r95_sys[1][np],
					shape_bys_r95_sys[2][np],
					shape_bys_r95_sys[3][np],
					shape_bys_r95_sys[4][np]
			      );
			printf("Bys(nosys)   %3.0f %10.5f %10.5f %10.5f %10.5f %10.5f\n",
					xReviewPoints[np],
					shape_bys_r95_nosys[0][np],
					shape_bys_r95_nosys[1][np],
					shape_bys_r95_nosys[2][np],
					shape_bys_r95_nosys[3][np],
					shape_bys_r95_nosys[4][np]
			      );
		}
		cout<<"\n\n"<<endl;
	}
	//----------end-----------------plot the projected exclusion limits 

	return 1;
}
