#ifndef  UNBINNEDBAYESIANLIMIT_H
#define  UNBINNEDBAYESIANLIMIT_H
#include <vector>
#include "BayesianLimitBase.h"
#include "Utilities.h"
#include "PdfRandom.h"
#include "CRandom.h"
using namespace std;
namespace lands{
class UnbinnedBayesianLimit
{
	public:
		UnbinnedBayesianLimit();
		~UnbinnedBayesianLimit();
		double P_vm_Given_brs(double r);		
		double Integral2(double a, double b, double fEpsilon=1.e-8);
		double RLimit(double alpha, double precision=1.e-4, double MinLikelihood=1.e-3, double integralPrecision=1.e-8, double rUpperBoud=1000.);
		void RunMCexps(int nexps, double alpha, double precision=1.e-4, double MinLikelihood=1.e-3, double integralPrecision=1.e-8, double rUpperBoud=1000.);
		vector<double> Get_vR(){return _vR;};
		double Get_R_mean();
		double Get_R_meanerr();
		void SetNpx(int npx){_npx=npx;};
		void SetRange(double xstart, double xstop){_xstart=xstart; _xstop=xstop;};
		void SetRdm(CRandom *rdm){_rdm=rdm;};

		//need a better name 
		double R_3A(double alpha, double precision=1.e-4, double MinLikelihood=1.e-3, double rUpperBoud=1000.);

		void SetPdfs(double (*fcn)(double *, double *), double *par); // set the pdf of signal distribution
		void SetPdfb(double (*fcn)(double *, double *), double *par);// set the pdf of bkg distribution
		void SetVmass(vector<double> v); //the unbinned input of mass
		void SetVmass(double *v, int nevts); //the unbinned input of mass
		void SetExpectedNsNb(double ns, double nb); // expected number of signal and bkg totally
		void SetDebug(bool d){debug=d;};
		double GetRLimit(){return _r95;};
	private:
		double (*_pdfs)(double *x, double *par); // normalized to unity
		double (*_pdfb)(double *x, double *par); // normalized to unity
		double *_pars; // parameters of signal pdf
		double *_parb; // parameters of bkg pdf

		double _ns;  // expected singal yield
		double _nb;  // expected bkg yield
		vector<double> _vmass; // mass list, unbinned
		bool debug;
		double _r95;

		CRandom *_rdm;
		PdfRandom *pdfRdmB, *pdfRdmS;
		double _xstart, _xstop;
		int _npx;
		vector<double> _vR;

};
};
#endif   /* ----- #ifndef UNBINNEDBAYESIANLIMIT_H  ----- */

