#ifndef BAYESIANLIMITBASE_H
#define BAYESIANLIMITBASE_H
#include <vector>
#include "CountingModel.h"
using namespace std;
namespace lands{
	// - --- - first step: simple counting experiment without systematic uncertainties
	// need interface for binned histograms as input 
	// check the precision of integral

	//============$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	//============$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	//		time consuming related parameters,  need to be configurable 
	//		and estimated when start a channel or combination calculation
	//const double rUpperBound = 1000; // fLL->Integral(0,rUpperBound) --> 1.0
	//const double precision = 1.e-4; //1.e-8
	//const double MinLikelihood = 1.e-3; //1.e-6
	//============$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
	//============$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

	const	double GreenBandLow = (1- 0.683)/2.; //1 sigma
	const	double GreenBandHigh = 1 - (1- 0.683)/2.; //1 sigma
	const	double YellowBandLow = (1- 0.955)/2.; //2 sigma
	const	double YellowBandHigh = 1 - (1- 0.955)/2.; //2 sigma

	double R_CL(CountingModel *cms,  double cl_alpha, double precision=1.e-30, 
			double MinLikelihood=1.e-6, double rUpperBoud=1000., bool debug=false); //calculate limit on sigma / sigma_SM , single channel or multiple channels(bins)
	double R_CL(double *par, double cl_alpha, int nbins, double precision=1.e-30, 
			double MinLikelihood=1.e-6, double rUpperBoud=1000., bool debug=false); //calculate limit on sigma / sigma_SM , single channel or multiple channels(bins)
	double R_CL(double nsig, double nbkg, double cl_alpha=0.95, double precision=1.e-30, double MinLikelihood=1.e-6, double rUpperBoud=1000., bool debug=false); //calculate limit on sigma / sigma_SM , single channel 
	double R_CL(double nsig, double nbkg, int ndat, double cl_alpha=0.95, double precision=1.e-30, double MinLikelihood=1.e-6, double rUpperBoud=1000., bool debug=false); //calculate limit on sigma / sigma_SM , single channel 
	double Likelihood_r(double *r, double *par, int npar); //calculate likelihood 
	double P_n0_Given_brs(double *r, double *par, int npar ); // the probability P(n0|b+rs)
	double P_n0_Given_brs(double *r, double *par); // the probability P(n0|b+rs)
	double R_CL_avr(double nsig, double nbkg, double cl_alpha); //single channel counting experiment, limit in average which is with mathmatical meaning

	double R_CL(vector<double> vs, vector<double> vb, double cl_alpha=0.95, double precision=1.e-30, double MinLikelihood=1.e-6, double rUpperBoud=1000., bool debug=false);
	double R_CL(vector<double> vs, vector<double> vb, vector<double> vd, double cl_alpha=0.95, double precision=1.e-30, double MinLikelihood=1.e-6, double rUpperBoud=1000., bool debug=false);
	double myIntegral(double (*fcn)(double *, double *, int), double a, double b, double *par, int npar, double steps=1000.);
};

#endif //#ifndef BAYESIANLIMITBASE_H
