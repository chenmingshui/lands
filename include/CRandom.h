#ifndef  CRANDOM_H
#define  CRANDOM_H
/*
// ported from ROOT TRandom3
 Random number generator class based on
   M. Matsumoto and T. Nishimura,
   Mersenne Twistor: A 623-diminsionally equidistributed
   uniform pseudorandom number generator
   ACM Transactions on Modeling and Computer Simulation,
   Vol. 8, No. 1, January 1998, pp 3--30.

 For more information see the Mersenne Twistor homepage
   http://www.math.keio.ac.jp/~matumoto/emt.html

 Advantage: large period 2**19937-1
            relativly fast
              (only two times slower than TRandom, but
               two times faster than TRandom2)
 Drawback:  a relative large internal state of 624 integers


 Aug.99 ROOT implementation based on CLHEP by P.Malzacher

 the original code contains the following copyright notice:
This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Library General Public
License as published by the Free Software Foundation; either
version 2 of the License, or (at your option) any later
version.
This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU Library General Public License for more details.
You should have received a copy of the GNU Library General
Public License along with this library; if not, write to the
Free Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
02111-1307  USA
Copyright (C) 1997 Makoto Matsumoto and Takuji Nishimura.
When you use this, send an email to: matumoto@math.keio.ac.jp
with an appropriate reference to your work.
*/
	typedef unsigned int UInt_t ;
	typedef double Double_t ;
	typedef int Int_t;
	typedef unsigned char UChar_t ; 
class CRandom{
	
	public:
	CRandom(UInt_t seed=4357);
	~CRandom();
	void SetSeed(UInt_t seed=0);
	double Rndm(int i=0);
	double Uniform(double xmin, double xmax);
	Int_t Poisson(double mean);
	double Gaus(double mean=0, double sigma=1);
	double Gamma(double mean=0);
	UInt_t GetSeed(){return fSeed;};
	
	private:
	Int_t fCount624;
	UInt_t fMt[624];
	UInt_t fSeed;

};

#endif   /* ----- #ifndef CRANDOM_H  ----- */

