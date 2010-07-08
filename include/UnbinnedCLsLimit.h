#ifndef  UNBINNEDCLSLIMIT_H
#define  UNBINNEDCLSLIMIT_H
#include <vector>
#include "Utilities.h"
#include "PdfRandom.h"
#include "CRandom.h"
#include "CLsLimit.h"
using namespace std;
namespace lands{
class UnbinnedCLsLimit
{
	public:
		UnbinnedCLsLimit();
		~UnbinnedCLsLimit();
		
		double lnQ(vector<double> vmass);
		double RLimit(double alpha=0.05, double epsilon=0.001, int nexps=100000);
		double RLimit(double alpha, double epsilon, double rmin, double rmax, int nexps=100000);
		double CLs(int nexps=100000);
		void SetPdfs(double (*fcn)(double *, double *), double *par); // set the pdf of signal distribution
		void SetPdfb(double (*fcn)(double *, double *), double *par);// set the pdf of bkg distribution
		void SetVmass(vector<double> v); //the unbinned input of mass
		void SetVmass(double *v, int nevts); //the unbinned input of mass
		void SetExpectedNsNb(double ns, double nb); // expected number of signal and bkg totally
		void SetDebug(bool d){debug=d;};
		void SetNpx(int npx){_npx=npx;};
		void SetRange(double xstart, double xstop){_xstart=xstart; _xstop=xstop;};
		void SetRdm(CRandom *rdm){_rdm=rdm;};
		vector<double> Get_vR(){return _vR;}; //
		vector<double> Get_vCLs(){return _vCLs;};//	
		vector<double> Get_m2logQ_b(){return _cls->Get_m2logQ_b();};
		vector<double> Get_m2logQ_sb(){return _cls->Get_m2logQ_sb();};
		double Get_m2logQ_d(){return -2*lnQ(_vmass);};
		bool CheckOk();
	private:
		double (*_pdfs)(double *x, double *par); // normalized to unity
		double (*_pdfb)(double *x, double *par); // normalized to unity
		double *_pars; // parameters of signal pdf
		double *_parb; // parameters of bkg pdf

		double _ns;  // expected singal yield
		double _nb;  // expected bkg yield
		vector<double> _vmass; // mass list, unbinned
		bool debug;

		PdfRandom *pdfRdmB, *pdfRdmS;
		double _xstart, _xstop;
		int _npx;

		CRandom *_rdm;

		vector<double> _vR;
		vector<double> _vCLs;

		CLsBase *_cls;
};
};
#endif   /* ----- #ifndef UNBINNEDBAYESIANLIMIT_H  ----- */

