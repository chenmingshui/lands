#include "CRandom.h"
#include <ctime>
#include <math.h>
#include <cmath>
#include <iostream>
// 
// ported from ROOT TRandom3
//

CRandom::CRandom(UInt_t seed){
	SetSeed(seed);
}
CRandom::~CRandom(){
}

Double_t CRandom::Rndm(Int_t)
{


	//  Machine independent random number generator.
	//  Produces uniformly-distributed floating points in ]0,1]
	//  Method: Mersenne Twistor

	UInt_t y;

	const Int_t  kM = 397;
	const Int_t  kN = 624;
	const UInt_t kTemperingMaskB =  0x9d2c5680;
	const UInt_t kTemperingMaskC =  0xefc60000;
	const UInt_t kUpperMask =       0x80000000;
	const UInt_t kLowerMask =       0x7fffffff;
	const UInt_t kMatrixA =         0x9908b0df;

	if (fCount624 >= kN) {
		register Int_t i;

		for (i=0; i < kN-kM; i++) {
			y = (fMt[i] & kUpperMask) | (fMt[i+1] & kLowerMask);
			fMt[i] = fMt[i+kM] ^ (y >> 1) ^ ((y & 0x1) ? kMatrixA : 0x0);
		}

		for (   ; i < kN-1    ; i++) {
			y = (fMt[i] & kUpperMask) | (fMt[i+1] & kLowerMask);
			fMt[i] = fMt[i+kM-kN] ^ (y >> 1) ^ ((y & 0x1) ? kMatrixA : 0x0);
		}

		y = (fMt[kN-1] & kUpperMask) | (fMt[0] & kLowerMask);
		fMt[kN-1] = fMt[kM-1] ^ (y >> 1) ^ ((y & 0x1) ? kMatrixA : 0x0);
		fCount624 = 0;
	}

	y = fMt[fCount624++];
	y ^=  (y >> 11);
	y ^= ((y << 7 ) & kTemperingMaskB );
	y ^= ((y << 15) & kTemperingMaskC );
	y ^=  (y >> 18);

	if (y) return ( (Double_t) y * 2.3283064365386963e-10); // * Power(2,-32)
	return Rndm();
}

Double_t CRandom::Uniform(Double_t x1, Double_t x2)
{
// returns a uniform deviate on the interval ]x1, x2].

   Double_t ans= Rndm();
   return x1 + (x2-x1)*ans;
}
void CRandom::SetSeed(UInt_t seed)
{
	//  Set the random generator sequence
	// if seed is 0 (default value) a TUUID is generated and used to fill
	// the first 8 integers of the seed array.
	// In this case the seed is guaranteed to be unique in space and time.
	// Use upgraded seeding procedure to fix a known problem when seeding with values
	// with many zero in the bit pattern (like 2**28).
	// see http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/MT2002/emt19937ar.html

	//   TRandom::SetSeed(seed);
	if( seed==0 ) {
		time_t curtime;      // Set 'random' seed number  if seed=0
		time(&curtime);      // Get current time in fSeed.
		fSeed = (UInt_t)curtime;
	} else {
		fSeed = seed;
	}

	fCount624 = 624;
	Int_t i,j;
	if (fSeed > 0) {
		fMt[0] = fSeed;
		j = 1;
	} else {
		std::cout<<"Warning, seed can't be set as =0.  Reset it to with random seed"<<std::endl;
		SetSeed(seed);
		/*
		   TUUID uid;
		   UChar_t uuid[16];
		   uid.GetUUID(uuid);
		   for (i=0;i<8;i++) {
		   fMt[i] = uuid[2*i]*256 +uuid[2*i+1];
		   if (i > 1) fMt[i] += fMt[0];
		   }
		   j = 8;
		 */
	}
	// use multipliers from  Knuth's "Art of Computer Programming" Vol. 2, 3rd Ed. p.106
	for(i=j; i<624; i++) {
		fMt[i] = (1812433253 * ( fMt[i-1]  ^ ( fMt[i-1] >> 30)) + i);
	}
}
Int_t CRandom::Poisson(Double_t mean)
{
	// Generates a random integer N according to a Poisson law.
	// Prob(N) = exp(-mean)*mean^N/Factorial(N)
	//
	// Use a different procedure according to the mean value.
	// The algorithm is the same used by CLHEP
	// For lower value (mean < 25) use the rejection method based on
	// the exponential
	// For higher values use a rejection method comparing with a Lorentzian
	// distribution, as suggested by several authors
	// This routine since is returning 32 bits integer will not work for values larger than 2*10**9
	// One should then use the Trandom::PoissonD for such large values
	//
	Int_t n;
	if (mean <= 0) return 0;
	if (mean < 25) {
		Double_t expmean = exp(-mean);
		Double_t pir = 1;
		n = -1;
		while(1) {
			n++;
			pir *= Rndm();
			if (pir <= expmean) break;
		}
		return n;
	}
	// for large value we use inversion method
	else if (mean < 1E9) {
		Double_t em, t, y;
		Double_t sq, alxm, g;
		Double_t pi = 3.14159265358979323846;

		sq = sqrt(2.0*mean);
		alxm = log(mean);
		g = mean*alxm - lgamma(mean + 1.0);

		do {
			do {
				y = tan(pi*Rndm());
				em = sq*y + mean;
			} while( em < 0.0 );

			em = floor(em);
			t = 0.9*(1.0 + y*y)* exp(em*alxm - lgamma(em + 1.0) - g);
		} while( Rndm() > t );

		return static_cast<Int_t> (em);

	}
	else {
		// use Gaussian approximation vor very large values
		// Gaus(0,1)

		Double_t arg = (0-1)/1;
		double gaus01 = exp(-0.5*arg*arg);

		n = Int_t(gaus01*sqrt(mean) + mean +0.5);
		return n;
	}
}
double CRandom::Gaus(double mean, double sigma){
//               
//  samples a random number from the standard Normal (Gaussian) Distribution 
//  with the given mean and sigma.                                                 
//  Uses the Acceptance-complement ratio from W. Hoermann and G. Derflinger 
//  This is one of the fastest existing method for generating normal random variables. 
//  It is a factor 2/3 faster than the polar (Box-Muller) method used in the previous 
//  version of TRandom::Gaus. The speed is comparable to the Ziggurat method (from Marsaglia)
//  implemented for example in GSL and available in the MathMore library. 
//                                                                           
//                                                                             
//  REFERENCE:  - W. Hoermann and G. Derflinger (1990):                       
//               The ACR Method for generating normal random variables,       
//               OR Spektrum 12 (1990), 181-185.                             
//                                                                           
//  Implementation taken from 
//   UNURAN (c) 2000  W. Hoermann & J. Leydold, Institut f. Statistik, WU Wien 
///////////////////////////////////////////////////////////////////////////////



   const Double_t kC1 = 1.448242853;
   const Double_t kC2 = 3.307147487;
   const Double_t kC3 = 1.46754004;
   const Double_t kD1 = 1.036467755;
   const Double_t kD2 = 5.295844968;
   const Double_t kD3 = 3.631288474;
   const Double_t kHm = 0.483941449;
   const Double_t kZm = 0.107981933;
   const Double_t kHp = 4.132731354;
   const Double_t kZp = 18.52161694;
   const Double_t kPhln = 0.4515827053;
   const Double_t kHm1 = 0.516058551;
   const Double_t kHp1 = 3.132731354;
   const Double_t kHzm = 0.375959516;
   const Double_t kHzmp = 0.591923442;
   /*zhm 0.967882898*/

   const Double_t kAs = 0.8853395638;
   const Double_t kBs = 0.2452635696;
   const Double_t kCs = 0.2770276848;
   const Double_t kB  = 0.5029324303;
   const Double_t kX0 = 0.4571828819;
   const Double_t kYm = 0.187308492 ;
   const Double_t kS  = 0.7270572718 ;
   const Double_t kT  = 0.03895759111;

   Double_t result;
   Double_t rn,x,y,z;


   do {
      y = Rndm();

      if (y>kHm1) {
         result = kHp*y-kHp1; break; }
  
      else if (y<kZm) {  
         rn = kZp*y-1;
         result = (rn>0) ? (1+rn) : (-1+rn);
         break;
      } 

      else if (y<kHm) {  
         rn = Rndm();
         rn = rn-1+rn;
         z = (rn>0) ? 2-rn : -2-rn;
         if ((kC1-y)*(kC3+fabs(z))<kC2) {
            result = z; break; }
         else {  
            x = rn*rn;
            if ((y+kD1)*(kD3+x)<kD2) {
               result = rn; break; }
            else if (kHzmp-y<exp(-(z*z+kPhln)/2)) {
               result = z; break; }
            else if (y+kHzm<exp(-(x+kPhln)/2)) {
               result = rn; break; }
         }
      }

      while (1) {
         x = Rndm();
         y = kYm * Rndm();
         z = kX0 - kS*x - y;
         if (z>0) 
            rn = 2+y/x;
         else {
            x = 1-x;
            y = kYm-y;
            rn = -(2+y/x);
         }
         if ((y-kAs+x)*(kCs+x)+kBs<0) {
            result = rn; break; }
         else if (y<x+kT)
            if (rn*rn<4*(kB-log(x))) {
               result = rn; break; }
      }
   } while(0);


   return mean + sigma * result;
}

double CRandom::Gamma(double a) {
/*  Joel Heinrich 25 March 2005
   Returns random deviate from gamma distribution.
   see D. Knuth "The Art of Computer Programming", vol 2, sec 3.4.1

   The p.d.f. of the gamma distribution is

              x^(a-1) e^(-x)
   f(x)dx  =  -------------- dx
                Gamma(a)

   where a, the order parameter, is the first argument of gamma_gen.

   The mean value and variance of the gamma distribution are both a.

   The user must supply a uniform random generator: gamma_gen's
   second argument is a pointer to a function that takes no arguments
   and returns a double, where the double is a uniform random deviate.

 */

  if (a >= 1.0e4) {
    //const double s=1/sqrt(9*a), x=1+s*(gauss_gen(ru)-s);
    const double s=1/sqrt(9*a), x=1+s*(Gaus()-s);
    return a*x*x*x; /* Wilson-Hilferty transformation */
  } else if (a > 4) {
    double x, y;
    const double aa=a-1 , scale = sqrt(2*a-1);
    do {
      x = scale * ( y = 1/tan(M_PI*Rndm()) );
    } while ( x <= -aa || Rndm() > (1+y*y)*exp(aa*log1p(x/aa)-x) ) ;
    return x+aa;
  } else {
    double dev=0;
    int ia = (int)a;
    const double fa = a-ia;
    if(ia>0) {
      double p = Rndm();
      while(--ia) p *= Rndm();
      dev = -log(p);
    }
    if(fa>0) {
      const double ra = 1/fa, rb = 1/(1-fa);
      /* Johnk's method for beta deviates employed here */
      for(;;) {
	const double u = pow(Rndm(),ra), w = u + pow(Rndm(),rb);
	if(w<=1) return dev - log(Rndm())*u/w;
      }
    }
    return dev;
  }
}
