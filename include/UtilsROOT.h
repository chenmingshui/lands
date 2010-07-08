#ifndef  UTILSROOT_H
#define  UTILSROOT_H
#include "TMath.h"
#include "Math/DistFunc.h"
#include "TTree.h"
#include <vector>

using namespace std;
/*
inline Double_t Significance(Double_t pvalue){
	// return sqrt(2.)*TMath::ErfInverse(1 - 2.*pvalue);
	// /afs/cern.ch/sw/lcg/app/releases/ROOT/5.26.00/src/root/math/mathcore/src/SpecFuncCephesInv.cxx
	return TMath::Abs(::ROOT::Math::normal_quantile(pvalue,1) ); 
}
*/
// FIXME make a few trees which save all important numbers, so we can easily remake plots 
// and parallel running .....,  combine them at the final step
void FillTree(TString sfile, vector<double> array);
void FillTree(TString sfile, vector<int> array);
void FillTree(TString sfile, double * array, int nsize=100000);
void FillTree(TString sfile, int* array, int nsize=100000);

#endif   /* ----- #ifndef UTILSROOT_INC  ----- */

