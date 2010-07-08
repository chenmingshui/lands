/*
 * =====================================================================================
 * 
 *       Filename:  PdfRandom.h
 * 
 *    Description:  generate random number following a pdf,  ported from root::TF1  
 * 
 * =====================================================================================
 */

#ifndef  PDFRANDOM_H
#define  PDFRANDOM_H
#include "CRandom.h"
class PdfRandom
{
//-- if one can provide the inverse funtion of the cumulative pdf, then we can throw random number faster.
//   will implement that later
	public:
		PdfRandom();
		~PdfRandom();
		void SetNpx(int npx=100);
		void Update();
		double GetRandom();
		double	GetRandom(double xmin, double xmax);
		void SetFunction(double (*fcn)(double *, double *), double *par, int npar=1);//normalized to unity
		double Integral(double a, double b, double epsilon=1.e-8);
		int BinarySearch(int n, const double  *array, double value); //independ on this class
		void SetRndGen(CRandom *rnd);
		void SetRange(double xmin, double xmax);
//		void SetParameters(double *par){fPar=par;};
		double* GetParameters(){ return fPar;};
		double GetParameter(int i){ return fPar[i];};
		void SetUgly(double ugly){_ugly=ugly;};
		double GetUgly(){return _ugly;};
	private:
		int fNpx;
		int fNpar;
		double fXmin;
		double fXmax;
		double* fAlpha;
		double* fBeta;
		double* fGamma;
		double* fIntegral;
		double  (*fFunction)(double *, double *);
		double* fPar;
		CRandom* fRdm;
		bool bHaveSetFunction;
		double _ugly;
	
};
#endif   /* ----- #ifndef PDFRANDOM_H  ----- */

