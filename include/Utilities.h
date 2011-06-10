#ifndef  UTILITIES_H
#define  UTILITIES_H
/*
 * =====================================================================================
 * 
 *       Filename:  Utilities.h
 * 
 *    Description:  
 * 
 *        Version:  1.0
 *        Created:  09/19/2009 10:36:50 PM CEST
 *       Revision:  none
 *       Compiler:  gcc
 * 
 *         Author:  Mingshui Chen (), Mingshui.Chen@cern.ch
 *        Company:  Univsity of Florida
 * 
 * =====================================================================================
 */

#include <algorithm>
#include <math.h>
#include <vector>

using std::vector;
using std::pair;
namespace lands{
//typedef bool bool;

template <typename Element, typename Index> void Sort(Index n, const Element* a, Index* index, bool down=true);

template<typename T> 
struct CompareDesc { 

   CompareDesc(T d) : fData(d) {}

   template<typename Index>
   bool operator()(Index i1, Index i2) { 
      return *(fData + i1) > *(fData + i2);
   }

   T fData;
};
template<typename T> 
struct CompareAsc { 

   CompareAsc(T d) : fData(d) {}

   template<typename Index>
   bool operator()(Index i1, Index i2) { 
      return *(fData + i1) < *(fData + i2);
   }

   T fData; 
};
template <typename Element, typename Index> void Sort(Index n, const Element* a, Index* index, bool down)
{
   // Sort the n elements of the  array a of generic templated type Element.
   // In output the array index of type Index contains the indices of the sorted array.
   // If down is false sort in increasing order (default is decreasing order).

   // NOTE that the array index must be created with a length >= n
   // before calling this function.
   // NOTE also that the size type for n must be the same type used for the index array
   // (templated type Index)

   for(Index i = 0; i < n; i++) { index[i] = i; }
   if ( down )
      std::sort(index, index + n, CompareDesc<const Element*>(a) );
   else
      std::sort(index, index + n, CompareAsc<const Element*>(a) );
}

double Poisson(double x, double par);  // function to calculate the probability of possion
double Poisson(double *x, double *par);
double Integral(double (*fcn)(double *, double *, int), double a, double b, double *par, int npar, double epsilon=1.e-30); //doing integration
double Integral(double (*fcn)(double *, double *), double a, double b, double *par, double epsilon=1.e-30);//doing integration
double FermiCurve(double *x, double *par);//to fit the turn on curve of P vs R
void erase(vector<double>& v, int pos);
void erase(vector<int>& v, int pos);
double LinearInterpolation(double x1, double y1, double x2, double y2, double y);
double LogLinearInterpolation(double x1, double y1, double x2, double y2, double y);
double LogLinearInterpolationErr(double x1, double y1, double ey1, double x2, double y2, double ey2, double y);
double FCInterpolation(double x1, double y1, double x2, double y2, double y);	
double IntegralSum(double (*fcn)(double *, double *, int), double a, double b, double *par, int npar);
int GetBandsByFermiCurveInterpolation(vector<double>  rn, vector<double> pn, double& _1SigmaLow, double& _1SigmaHigh, double& _2SigmaLow, double& _2SigmaHigh);
int GetBandsByLinearInterpolation(vector<double>  rn, vector<double> pn, double& _1SigmaLow, double& _1SigmaHigh, double& _2SigmaLow, double& _2SigmaHigh);
int GetBandsByNoInterpolation(vector<double>  rn, vector<double> pn, double& _1SigmaLow, double& _1SigmaHigh, double& _2SigmaLow, double& _2SigmaHigh, double &median);
int GetBandsByFeldmanCousins(vector<double>  rn, vector<double> pn, double& _1SigmaLow, double& _1SigmaHigh, double& _2SigmaLow, double& _2SigmaHigh);
int GetBands(vector<double> & rn, vector<double>& pn, double& _1SigmaLow, double& _1SigmaHigh, double& _2SigmaLow, double& _2SigmaHigh);
int GetBands(vector<double> & vx, double& _1SigmaLow, double& _1SigmaHigh, double& _2SigmaLow, double& _2SigmaHigh);
int GetBands(double *dx, int nexps, double& _1SigmaLow, double& _1SigmaHigh, double& _2SigmaLow, double& _2SigmaHigh);
double GetBandByFermiCurveInterpolation(vector<double>  rn, vector<double> pn,  double sigma);
double GetBandByLinearInterpolation(vector<double>  rn, vector<double> pn,  double sigma);
void SortAndCumulative(double *Q_b, int nexps, vector<double> & qn, vector<double>& pn, bool down=0); // by default, it's sorted by increased order
void SortAndCumulative(vector<double> Q_b, vector<double> & qn, vector<double>& pn, bool down=0);
double GetMeanOfSortedXwithProb(vector<double> vx, vector<double> vp );

double TruncatedGaussianPdf(double *x, double *par);


	// n*(n-1)*(n-2)*.......*2*1
	long nto1(int i);  // this is Gamma Function
	// single channel counting experiment  CLs,  identical result with Bayesian approach, mathmaticlly
	double CLs_Analytics(double s, double b, int d);  // with data = integer
	double CLs_Analytics(double s, double b);		// with data = b, not integer  // using continuous Poisson 
	double CLs_Analytics(double s, double b, double d);    // with data not integer   // using continuous Poisson
	vector< pair<double, double> > m2lnQ_b_Analytics(double s, double b);
	vector< pair<double, double> > m2lnQ_sb_Analytics(double s, double b);
	double m2lnQ(double s,double b, double d);

	double fP_n0_Given_brs(double *r, double *par,int npar);

	void gausslaguerre(double x[],double lw[],int n,double alpha);

	namespace Cephes { 
		// special functions taken from Cephes library 
		//  see:  http://www.netlib.org/cephes
		// 
		// Copyright 1985, 1987, 2000 by Stephen L. Moshier
		// 
		//  granted permission from the author to be used in MathCore
		//  

		/* routines for efficient polynomial evaluation*/
		double Polynomialeval(double x, double* a, unsigned int N);
		double Polynomial1eval(double x, double* a, unsigned int N);


		//---
		/* the machine roundoff error */
#define kMACHEP  1.11022302462515654042363166809e-16

		/* largest argument for TMath::Exp() */
#define kMAXLOG  709.782712893383973096206318587

		/* smallest argument for TMath::Exp() without underflow */
#define kMINLOG  -708.396418532264078748994506896

		/* the maximal number that pow(x,x-0.5) has no overflow */
		/* we use a (very) conservative portable bound          */
#define kMAXSTIR  108.116855767857671821730036754

#define kMAXLGM 2.556348e305

		/* normal quantile */ 
		double ndtri (double y); 

	} // end namespace Cephes

double normal_quantile(double z, double sigma) ;
inline double Significance(double pvalue){
	// return sqrt(2.)*TMath::ErfInverse(1 - 2.*pvalue);
	// /afs/cern.ch/sw/lcg/app/releases/ROOT/5.26.00/src/root/math/mathcore/src/SpecFuncCephesInv.cxx
	//return TMath::Abs(::ROOT::Math::normal_quantile(pvalue,1) ); 
	return fabs(normal_quantile(pvalue,1) ); 
}
double InverseCDF(
		vector<double> v,
		double alpha,
		double pvalue, 
		double sigmaVariation, 
		double& inverseWithVariation);
};
#endif   /* ----- #ifndef UTILITIES_H  ----- */

