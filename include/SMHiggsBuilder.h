#ifndef SMHiggsBuilder_H
#define SMHiggsBuilder_H
#include <iostream>
#include <cmath>
#include <string>
#include <cstdlib>
#include <fstream>

#include <Math/Interpolator.h>
#include "TROOT.h"
#include "TF1.h"
#include "TGraph.h"
#include "TSpline.h"
#include "TString.h"

using namespace std;
namespace lands{
	class SMHiggsBuilder
	{
		public:
			SMHiggsBuilder();
			~SMHiggsBuilder();
			//void init();
			void readSMBr(TString sf);
			double br(int decayMode,double mH);

		private:
			double mass_BR[217];
			double BR[26][217];
			int N_BR;
			std::string FileLoc;
			vector<ROOT::Math::Interpolator*> vsl;  
	};
}
#endif   /* ----- #ifndef SMHiggsBuilder_h----- */
