#include "SMHiggsBuilder.h"
namespace lands{
	SMHiggsBuilder::SMHiggsBuilder(){
		vsl.clear();
	}
	SMHiggsBuilder::~SMHiggsBuilder(){
		for(int i=0; i<vsl.size(); i++) {if(vsl[i]) delete vsl[i];}
		vsl.clear();
	}
	//void SMHiggsBuilder::init(){readSMBr();}
	void SMHiggsBuilder::readSMBr(TString sf){ 
		N_BR = 217;
		ifstream file;

		// Read Widths into memory
		if(sf=="")
		FileLoc = "../data/HiggsBR_7TeV_Official.txt"; //directory of input file
		else FileLoc = sf;
		//cout<<" Reading SM Branch ratio from "<<FileLoc.c_str()<<endl;
		const char* BranchRatioFileLoc = FileLoc.c_str(); 
		file.open(BranchRatioFileLoc);
		if(!file) {cout<<" SM Br file "<<BranchRatioFileLoc<<" doesn't exist !"<<endl; exit(1);}
		for(int k = 0; k < N_BR; k++){

			file >> mass_BR[k] >> BR[0][k] >> BR[1][k] >> BR[2][k] >> BR[3][k] >> BR[4][k] >> BR[5][k] >> BR[6][k] >> BR[7][k] >> BR[8][k] >> BR[9][k]
				>> BR[10][k] >> BR[11][k] >> BR[12][k] >> BR[13][k] >> BR[14][k] >> BR[15][k] >> BR[16][k] >> BR[17][k] >> BR[18][k] >> BR[19][k] >> BR[20][k]
				>> BR[21][k] >> BR[22][k] >> BR[23][k] >> BR[24][k] >> BR[25][k];


		}
		file.close();


		/***********************IDs************************/
		/*                       Total = 0                */
		/*                       H->bb = 1                */
		/*                   H->tautau = 2                */
		/*                     H->mumu = 3                */
		/*                       H->ss = 4                */
		/*                       H->cc = 5                */
		/*                       H->tt = 6                */
		/*                       H->gg = 7                */
		/*                   H->gamgam = 8                */
		/*                     H->gamZ = 9                */
		/*                       H->WW = 10               */
		/*                       H->ZZ = 11               */
		/*                       H->4e = 12               */
		/*                    H->2e2mu = 13               */
		/*              H->4lep (e,mu) = 14               */
		/*          H->4lep (e,mu,tau) = 15               */
		/*                H->e+nu e-nu = 16               */
		/*               H->e+nu mu-nu = 17               */
		/*    H->2l2nu(l=e,mu)(nu=any) = 18               */
		/* H->2l2nu(l=e,mu,tau)(nu=any) = 19              */  
		/*    H->2l2q (l=e,mu)(q=udcsb) = 20              */
		/* H->2l2q(l=e,mu,tau)(q=udcsb) = 21              */
		/* H->l+nu qq(*) (l=e,mu)(q=udcsb) = 22           */
		/*  H->2nu2q (nu=any)(q=udcsb) = 23               */
		/*            H->4q (q=udcsb) = 24                */
		/*      H->4f (f=any fermion) = 25                */
		/**************************************************/

		for(int i=0; i<26; i++){
			ROOT::Math::Interpolator* interp = new ROOT::Math::Interpolator(N_BR, ROOT::Math::Interpolation::kCSPLINE);
			interp->SetData(N_BR, mass_BR, BR[i]);
			vsl.push_back(interp);
		}

	}
	double SMHiggsBuilder::br(int dm, double mH){
		return vsl[dm]->Eval(mH);
	}
}
