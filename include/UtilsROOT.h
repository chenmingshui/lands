#ifndef  UTILSROOT_H
#define  UTILSROOT_H
#include "TMath.h"
#include "Math/DistFunc.h"
#include "TTree.h"
#include "CountingModel.h"
#include "TString.h"
#include <vector>
#include <string>
#include "TH1F.h"
#include "TROOT.h"
#include "TGraphErrors.h"

using namespace std;
using namespace lands;
/*
inline Double_t Significance(Double_t pvalue){
	// return sqrt(2.)*TMath::ErfInverse(1 - 2.*pvalue);
	// /afs/cern.ch/sw/lcg/app/releases/ROOT/5.26.00/src/root/math/mathcore/src/SpecFuncCephesInv.cxx
	return TMath::Abs(::ROOT::Math::normal_quantile(pvalue,1) ); 
}
*/
// FIXME make a few trees which save all important numbers, so we can easily remake plots 
// and parallel running .....,  combine them at the final step
void SaveResults(TString sfile, double mH, double limit, double limitErr, double significance, double pvalue, double rm2s, double rm1s, double rmedian, double rmean, double rp1s, double rp2s);
void FillTree(TString sfile, vector<double> array, TString treeName="");
void FillTree(TString sfile, vector<int> array);
void FillTree(TString sfile, double * array, int nsize=100000);
void FillTree(TString sfile, int* array, int nsize=100000);
void FillTree(TString sfile, double d1, double d2, vector<double> array1,  vector<double> array2, TString d1Name="d1", TString d2Name="d2", TString array1Name="T1", TString array2Name="T2", TString option = "RECREATE");
void FillTree2(TString sfile, double d1, double d2, vector<double> array1,  vector<double> array2, TString d1Name="d1", TString array1Name="T1", TString array2Name="T2", TString option = "RECREATE");

bool isWordInMap(TString s, std::map<TString, vector<TString> > tMap);
void StringStrip( std::string & str ) ;
void StringSplit( std::vector < std::string > & splitValues, 
		const std::string & str,
		const std::string & delim ); 
TString ReadFile(const char*fileName);
bool ConfigureModel(CountingModel *cms, double mass,  TString ifileContentStripped, int debug=0);
bool ConfigureModel(CountingModel *cms, double mass, const char* fileName, int debug=0);
bool CheckIfDoingShapeAnalysis(CountingModel* cms, double mass, TString ifileContentStripped, int debug=0);
vector<TString> SplitIntoLines(TString ifileContentStripped, bool debug=false);
TString GetWordFromLine(TString line, int index, string delim = " ");

TTree * LoadTreeBonly(TString filename, TString & treeName);
TH1F* GetHisto(string filename, string histoname);
TObject* GetTObject(string filename, string objname);
bool ConfigureShapeModel(CountingModel *cms, double mass, TString ifileContentStripped, vector< vector<string> > parametricShapeLines,  vector< vector<string> > uncerlinesAffectingShapes,  int debug=0);
RooAbsData* GetRooAbsData(string c, string p, vector< vector<string> > lines, double mass=0);
RooAbsPdf* GetPdf(string c, string p, vector< vector<string> > lines, double mass=0);
RooAbsArg* GetExtraNorm(string c, string p, vector< vector<string> > lines, double mass=0);

void ReadLimitVsCLsFromFile(TGraphErrors*tge, TFile*f, int debug=0); 
bool GetCLs(double qdata, TTree* tsb, TTree*tb,  double &cls, double &err, int debug=0);
bool GetM2lnQ(TTree* tsb, TTree*tb, vector<double> &vclsb, vector<double>&vclb, int debug=0);
bool GetPValue(vector<double> vclsb, double qdata, double &ret, double &err, int debug=0);
void ReadM2lnQGridFromFile(TString filename, std::map<double, TTree*>&gridCLsb, std::map<double, TTree*>&gridCLb, int _debug=0); 
vector<double> GetVectorFrom(TTree* tree, TString brName);
void ReadM2lnQGridFromFile(TString filename, std::map<double, double>&gridQdata, int _debug=0); 
#endif   /* ----- #ifndef UTILSROOT_INC  ----- */

