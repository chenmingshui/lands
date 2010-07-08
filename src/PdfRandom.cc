// @(#)root/hist:$Id: TF1.h 29895 2009-08-25 09:35:06Z brun $
// Author: Rene Brun   18/08/95
/*
 * =====================================================================================
 *
 *       Filename:  PdfRandom.cc   ported from ROOT TF1
 * =====================================================================================
 */
#include "PdfRandom.h"
#include <math.h>
#include <iostream>
PdfRandom::PdfRandom(){
   //fXmin      = 0;
   //fXmax      = 0;
   fNpx       = 100;
   fIntegral  = 0;
   fAlpha     = 0;
   fBeta      = 0;
   fGamma     = 0;
   fPar	      = 0;
	bHaveSetFunction=false;
}
PdfRandom::~PdfRandom()
{
   if (fIntegral)  delete [] fIntegral;
   if (fAlpha)     delete [] fAlpha;
   if (fBeta)      delete [] fBeta;
   if (fGamma)     delete [] fGamma;
   if (fPar)	   delete [] fPar;
//	fPar=NULL;
	fRdm=NULL;
}
double PdfRandom::GetRandom()
{
   // Return a random number following this function shape
   //
   //   The distribution contained in the function fname (TF1) is integrated
   //   over the channel contents.
   //   It is normalized to 1.
   //   For each bin the integral is approximated by a parabola.
   //   The parabola coefficients are stored as non persistent data members
   //   Getting one random number implies:
   //     - Generating a random number between 0 and 1 (say r1)
   //     - Look in which bin in the normalized integral r1 corresponds to
   //     - Evaluate the parabolic curve in the selected bin to find
   //       the corresponding X value.
   //   if the ratio fXmax/fXmin > fNpx the integral is tabulated in log scale in x
   //   The parabolic approximation is very good as soon as the number
   //   of bins is greater than 50.

   //  Check if integral array must be build
	if(!bHaveSetFunction) return 0;

   if (fIntegral == 0) {
      fIntegral = new double[fNpx+1];
      fAlpha    = new double[fNpx+1];
      fBeta     = new double[fNpx];
      fGamma    = new double[fNpx];
      fIntegral[0] = 0;
      fAlpha[fNpx] = 0;
      double integ;
      int intNegative = 0;
      int i;
      bool logbin = false;
      double dx;
      double xmin = fXmin;
      double xmax = fXmax;
      if (xmin > 0 && xmax/xmin> fNpx) {
         logbin =  true;
         fAlpha[fNpx] = 1;
         xmin = log10(fXmin);
         xmax = log10(fXmax);
      }
      dx = (xmax-xmin)/fNpx;
         
      double *xx = new double[fNpx+1];
      for (i=0;i<fNpx;i++) {
            xx[i] = xmin +i*dx;
      }
      xx[fNpx] = xmax;
      for (i=0;i<fNpx;i++) {
         if (logbin) {
            integ = Integral(pow(10,xx[i]), pow(10,xx[i+1]));
         } else {
            integ = Integral(xx[i],xx[i+1]);
         }
         if (integ < 0) {intNegative++; integ = -integ;}
         fIntegral[i+1] = fIntegral[i] + integ;
      }
      if (intNegative > 0) {
         std::cout<<"Warning::GetRandom function: has  negative values: abs assumed"<<std::endl;
      }
      if (fIntegral[fNpx] == 0) {
         std::cout<<"Error::GetRandom  Integral of function is zero"<<std::endl;
         return 0;
      }
      double total = fIntegral[fNpx];
      for (i=1;i<=fNpx;i++) {  // normalize integral to 1
         fIntegral[i] /= total;
      }
      //the integral r for each bin is approximated by a parabola
      //  x = alpha + beta*r +gamma*r**2
      // compute the coefficients alpha, beta, gamma for each bin
      double x0,r1,r2,r3;
      for (i=0;i<fNpx;i++) {
         x0 = xx[i];
         r2 = fIntegral[i+1] - fIntegral[i];
         if (logbin) r1 = Integral(pow(10,x0),pow(10,x0+0.5*dx))/total;
         else        r1 = Integral(x0,x0+0.5*dx)/total;
         r3 = 2*r2 - 4*r1;
         if (fabs(r3) > 1e-8) fGamma[i] = r3/(dx*dx);
         else           fGamma[i] = 0;
         fBeta[i]  = r2/dx - fGamma[i]*dx;
         fAlpha[i] = x0;
         fGamma[i] *= 2;
      }
      delete [] xx;
   }

   // return random number
   double r  = fRdm->Rndm();
   int bin  = BinarySearch(fNpx,fIntegral,r);
   double rr = r - fIntegral[bin];

   double yy;
   if(fGamma[bin] != 0)
      yy = (-fBeta[bin] + sqrt(fBeta[bin]*fBeta[bin]+2*fGamma[bin]*rr))/fGamma[bin];
   else
      yy = rr/fBeta[bin];
   double x = fAlpha[bin] + yy;
   if (fAlpha[fNpx] > 0) return pow(10,x);
   return x;
}


//______________________________________________________________________________
double PdfRandom::GetRandom(double xmin, double xmax)
{
   // Return a random number following this function shape in [xmin,xmax]
   //
   //   The distribution contained in the function fname (TF1) is integrated
   //   over the channel contents.
   //   It is normalized to 1.
   //   For each bin the integral is approximated by a parabola.
   //   The parabola coefficients are stored as non persistent data members
   //   Getting one random number implies:
   //     - Generating a random number between 0 and 1 (say r1)
   //     - Look in which bin in the normalized integral r1 corresponds to
   //     - Evaluate the parabolic curve in the selected bin to find
   //       the corresponding X value.
   //   The parabolic approximation is very good as soon as the number
   //   of bins is greater than 50.
   //
   //  IMPORTANT NOTE
   //  The integral of the function is computed at fNpx points. If the function
   //  has sharp peaks, you should increase the number of points (SetNpx)
   //  such that the peak is correctly tabulated at several points.

   //  Check if integral array must be build
   if (fIntegral == 0) {
      double dx = (fXmax-fXmin)/fNpx;
      fIntegral = new double[fNpx+1];
      fAlpha    = new double[fNpx];
      fBeta     = new double[fNpx];
      fGamma    = new double[fNpx];
      fIntegral[0] = 0;
      double integ;
      int intNegative = 0;
      int i;
      for (i=0;i<fNpx;i++) {
         integ = Integral(double(fXmin+i*dx), double(fXmin+i*dx+dx));
         if (integ < 0) {intNegative++; integ = -integ;}
         fIntegral[i+1] = fIntegral[i] + integ;
      }
      if (intNegative > 0) {
         std::cout<<"Warning::GetRandom function: has  negative values: abs assumed"<<std::endl;
      }
      if (fIntegral[fNpx] == 0) {
         std::cout<<"Error::GetRandom  Integral of function is zero"<<std::endl;
         return 0;
      }
      double total = fIntegral[fNpx];
      for (i=1;i<=fNpx;i++) {  // normalize integral to 1
         fIntegral[i] /= total;
      }
      //the integral r for each bin is approximated by a parabola
      //  x = alpha + beta*r +gamma*r**2
      // compute the coefficients alpha, beta, gamma for each bin
      double x0,r1,r2,r3;
      for (i=0;i<fNpx;i++) {
         x0 = fXmin+i*dx;
         r2 = fIntegral[i+1] - fIntegral[i];
         r1 = Integral(x0,x0+0.5*dx)/total;
         r3 = 2*r2 - 4*r1;
         if (fabs(r3) > 1e-8) fGamma[i] = r3/(dx*dx);
         else           fGamma[i] = 0;
         fBeta[i]  = r2/dx - fGamma[i]*dx;
         fAlpha[i] = x0;
         fGamma[i] *= 2;
      }
   }

   // return random number
   double dx   = (fXmax-fXmin)/fNpx;
   int nbinmin = (int)((xmin-fXmin)/dx);
   int nbinmax = (int)((xmax-fXmin)/dx)+2;
   if(nbinmax>fNpx) nbinmax=fNpx;

   double pmin=fIntegral[nbinmin];
   double pmax=fIntegral[nbinmax];

   double r,x,xx,rr;
   do {
      r  = fRdm->Uniform(pmin,pmax);

      int bin  = BinarySearch(fNpx,fIntegral,r);
      rr = r - fIntegral[bin];

      if(fGamma[bin] != 0)
         xx = (-fBeta[bin] + sqrt(fBeta[bin]*fBeta[bin]+2*fGamma[bin]*rr))/fGamma[bin];
      else
         xx = rr/fBeta[bin];
      x = fAlpha[bin] + xx;
   } while(x<xmin || x>xmax);
   return x;
}
void PdfRandom::Update()
{
   // Called by functions such as SetRange, SetNpx, SetParameters
   // to force the deletion of the associated histogram or Integral
   if (fIntegral) {
      delete [] fIntegral; fIntegral = 0;
      delete [] fAlpha;    fAlpha    = 0;
      delete [] fBeta;     fBeta     = 0;
      delete [] fGamma;    fGamma    = 0;
   }
}
void PdfRandom::SetNpx(int npx)
{
   // Set the number of points used to draw the function
   //
   // The default number of points along x is 100 for 1-d functions and 30 for 2-d/3-d functions
   // You can increase this value to get a better resolution when drawing
   // pictures with sharp peaks or to get a better result when using TF1::GetRandom
   // the minimum number of points is 4, the maximum is 100000 for 1-d and 10000 for 2-d/3-d functions

   if (npx < 4) {
      std::cout<<"Warning::SetNpx  Number of points must be >4 && < 100000, fNpx set to 4"<<std::endl;
      fNpx = 4;
   } else if(npx > 100000) {
      std::cout<<"Warning::SetNpx  Number of points must be >4 && < 100000, fNpx set to 100000"<<std::endl;
      fNpx = 100000;
   } else {
      fNpx = npx;
   }
   Update();
}
void PdfRandom::SetFunction(double (*fcn)(double *, double *), double *par, int npar){//normalized to unity
	fFunction=fcn;
	fNpar=npar;
	fPar=new double[fNpar];//FIXME maximum number of parameters
	for(int i=0; i<fNpar; i++)
		fPar[i]=par[i];
	bHaveSetFunction=true;
}
double PdfRandom::Integral(double a, double b, double epsilon)
{
	double fEpsilon = epsilon; // configurable  
	bool  fgAbsValue = true;

	const double kHF = 0.5;
	const double kCST = 5./1000;

	double x[12] = { 0.96028985649753623,  0.79666647741362674,
		0.52553240991632899,  0.18343464249564980,
		0.98940093499164993,  0.94457502307323258,
		0.86563120238783174,  0.75540440835500303,
		0.61787624440264375,  0.45801677765722739,
		0.28160355077925891,  0.09501250983763744};

	double w[12] = { 0.10122853629037626,  0.22238103445337447,
		0.31370664587788729,  0.36268378337836198,
		0.02715245941175409,  0.06225352393864789,
		0.09515851168249278,  0.12462897125553387,
		0.14959598881657673,  0.16915651939500254,
		0.18260341504492359,  0.18945061045506850};

	double h, aconst, bb, aa, c1, c2, u, s8, s16, f1, f2;
	double xx[1];
	int i;

	/*   if ( fFunction == 0 )
		 {
		 MATH_ERROR_MSG("ROOT::Math::GausIntegratorOneDim", "A function must be set first!");
		 return 0.0;
		 }
	 */
	h = 0;
	if (b == a) return h;
	if(b>a)aconst = kCST/(b-a);
	else aconst = - kCST/(b-a);

	bb = a;
CASE1:
	aa = bb;
	bb = b;
CASE2:
	c1 = kHF*(bb+aa);
	c2 = kHF*(bb-aa);
	s8 = 0;
	for (i=0;i<4;i++) {
		u     = c2*x[i];
		xx[0] = c1+u;
		//      f1    = (*fFunction)(xx);
		// --Mingshui
		//cout<<xx[0]<<" "<<par[0]<<" "<<par[1]<<endl;
		f1 = (*fFunction)(xx, fPar);


		if (fgAbsValue) {if(f1<0)f1 = -f1;}
		xx[0] = c1-u;
		//f2    = (*fFunction) (xx);
		// --Mingshui
		f2 = (*fFunction)(xx, fPar);

		if (fgAbsValue) {if(f2<0)f2 = -f2;}
		s8   += w[i]*(f1 + f2);
	}
	s16 = 0;
	for (i=4;i<12;i++) {
		u     = c2*x[i];
		xx[0] = c1+u;
		//f1    = (*fFunction) (xx);
		// --Mingshui
		f1 = (*fFunction)(xx, fPar);

		if (fgAbsValue) {if(f1<0)f1 = -f1;}
		xx[0] = c1-u;
		//f2    = (*fFunction) (xx);
		// --Mingshui
		f2 = (*fFunction)(xx, fPar);

		if (fgAbsValue) {if(f2<0)f2 = -f2;}
		s16  += w[i]*(f1 + f2);
	}
	s16 = c2*s16;
	double s16_tmp = s16;
	double s16_c2s8_tmp = s16-c2*s8;
	double c2_tmp=c2;


	if(s16_c2s8_tmp<0) s16_c2s8_tmp= -s16_c2s8_tmp;
	if(s16_tmp<0) s16_tmp=-s16_tmp;
	if(c2_tmp<0) c2_tmp=-c2_tmp;

	if (s16_c2s8_tmp <= fEpsilon*(1. + s16_tmp)) {
		h += s16;
		if(bb != b) goto CASE1;
	} else {
		bb = c1;
		if(1. + aconst*c2_tmp != 1) goto CASE2;
		h = s8;  //this is a crude approximation (cernlib function returned 0 !)
	}

	//  fUsedOnce = true;
	//  fLastResult = h;
	//  fLastError = std::abs(s16-c2*s8);

	return h;
}
int PdfRandom::BinarySearch(int n, const double  *array, double value)
{
	// Binary search in an array of n values to locate value.
	//
	// Array is supposed  to be sorted prior to this call.
	// If match is found, function returns position of element.
	// If no match found, function gives nearest element smaller than value.

	const double* pind;
	pind = std::lower_bound(array, array + n, value);
	if ( (pind != array + n) && (*pind == value) )
		return (pind - array);
	else
		return ( pind - array - 1);
}
void PdfRandom::SetRndGen(CRandom *rnd){
	fRdm=rnd;
}
void PdfRandom::SetRange(double xmin, double xmax)
{
	// Initialize the upper and lower bounds to draw the function.
	//
	// The function range is also used in an histogram fit operation
	// when the option "R" is specified.

	fXmin = xmin;
	fXmax = xmax;
	Update();
}
