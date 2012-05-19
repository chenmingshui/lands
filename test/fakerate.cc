#include "TMinuit.h"
#include <iostream>
#include <iomanip>
#include <math.h>
#include <ctime> // upto second
#include <utility>
#include "TMath.h"
#include "TPaveText.h"
#include <vector>
#include "PlotUtilities.h"
#include "UtilsROOT.h"
#include "Utilities.h"
using std::cout;
using std::endl;

	// Andrey Korytov's  hand-written model description  	
	// ****** user provide following information ***************
/*	int N_LL = 494.8; // Loose-Loose control region counts
	int N_TL = 29.95; // Tight-Loose control region counts
	float eff = 0.021; // fake rate  (must  be < 1)
	float eff_err = 0.30; // 20% uncertainty on the measured fake rate, assuming LogNormal pdf
*/
	int N_LL = 19; // Loose-Loose control region counts
	int N_TL = 8; // Tight-Loose control region counts
	float eff = 0.16; // fake rate  (must  be < 1)
	float eff_err = 0.30; // 20% uncertainty on the measured fake rate, assuming LogNormal pdf
	
	// if need to discriminate TL and LT, eg. emu and mue, second lepton as fake
	int N_eLmT= 20; // Tight-Loose (em) control region counts
	int N_mLeT= 20; // Loose-Tight (em) control region counts
	float effe = 0.1;
	float effm = 0.2;
	float effe_err = 0.3; // fraction error 
	float effm_err = 0.4;
	
	// ****** end of input information ***************
void Function_A(Int_t &npar, Double_t *gin, Double_t &f,  Double_t *par, Int_t iflag){	

	//float r2 = N_LL/(1-eff)/(1-eff); // i.e. trueFF,  initial estimation of the true FakeFake events, will be treated as flat pdf,  >=0

	// POI , i.e. R0 ,  the estimated total fake bkgs in signal region  
	// trueFP, i.e. r1,   true FakePrompt events
	
	float R0 = par[0];    // POI
	float r2 = par[1];    // flat
	float geff = par[2];  // underline nuisance of eff_err,   unit gaussian


	// L = R2^NLL * e^(-R2) / NLL!  *  R1^NTL * e^(-R1) / NTL!  
	// need also a term for Gaussian constraint on eff_err, i.e.  1/sqrt(2*PI)*exp(-x^2/2)

	float effcor = eff * pow(1+eff_err,geff);
	float R2 = r2*(1-effcor)*(1-effcor);
	float r1 = ( R0 - r2*effcor*effcor )/effcor;
	//float r1 = R0; // ( R0 - r2*effcor*effcor )/effcor;
	if(r1 < 0 ) {
		f=9e5; return;
	}


	float R1 = 2*r2*effcor*(1-effcor) + r1*(1-effcor) ;

	// lnL = NLL * R2 - R2 - ln(NLL!) + NTL * R1 - R1 - ln(NTL!)
	// ln(NLL!) = TMath::LnGamma(NLL+1);
	
	float lnL = N_LL * log(R2) - R2 - TMath::LnGamma(N_LL+1) + N_TL*log(R1) - R1 - TMath::LnGamma(N_TL+1);	
	lnL -= geff*geff/2.;

	f =  -2 * lnL;

//	par[3] = r1*effcor + r2*effcor*effcor;
}

void Function_B(Int_t &npar, Double_t *gin, Double_t &f,  Double_t *par, Int_t iflag){	

	//float r2 = N_LL/(1-eff)/(1-eff); // i.e. trueFF,  initial estimation of the true FakeFake events, will be treated as flat pdf,  >=0

	// POI , i.e. R0 ,  the estimated total fake bkgs in signal region  
	// trueFP, i.e. r1,   true FakePrompt events
	
	float R0 = par[0];    // POI
	float r2 = par[1];    // flat  >=0
	float reLmT = par[2];    // flat  >=0
	float geffe  = par[3];  // underline nuisance of eff_err,   unit gaussian
	float geffm  = par[4];  // underline nuisance of eff_err,   unit gaussian


	// L = R2^NLL * e^(-R2) / NLL!  *  R1^NTL * e^(-R1) / NTL!  
	// need also a term for Gaussian constraint on eff_err, i.e.  1/sqrt(2*PI)*exp(-x^2/2)

	float effecor = effe * pow(1+effe_err,geffe);
	float effmcor = effm * pow(1+effm_err,geffm);
	float R2 = r2*(1-effecor)*(1-effmcor);
	float rmLeT = ( R0 - r2*effecor*effmcor - reLmT*effecor)/effmcor;
	//float r1 = R0; // ( R0 - r2*effcor*effcor )/effcor;
	if(rmLeT < 0 ) {
		f=9e5; return;
	}


	float ReLmT = r2*effmcor*(1-effecor) + reLmT*(1-effecor) ;
	float RmLeT = r2*effecor*(1-effmcor) + rmLeT*(1-effmcor) ;

	// lnL = NLL * R2 - R2 - ln(NLL!) + NTL * R1 - R1 - ln(NTL!)
	// ln(NLL!) = TMath::LnGamma(NLL+1);
	
	float lnL = N_LL * log(R2) - R2 - TMath::LnGamma(N_LL+1) ;
	lnL += N_eLmT*log(ReLmT) - ReLmT - TMath::LnGamma(N_eLmT+1);	
	lnL += N_mLeT*log(RmLeT) - RmLeT - TMath::LnGamma(N_mLeT+1);	
	lnL -= geffe*geffe/2.;
	lnL -= geffm*geffm/2.;

	f =  -2 * lnL;

}
double interpolaton(TGraph *g, double y, double xbest, TString side = "left");
int main(int argc, const char*argv[]){

	int NumberOfKindsOfFakeObjects = 1;
	if(argc>2) {
		NumberOfKindsOfFakeObjects = std::atoi(argv[1]);
		if(NumberOfKindsOfFakeObjects==1){
			if(argc!=6) {
				cout<<" Usage:  ./fakerate.exe 1 N_LL N_TL FakeRate FractionError"<<endl; exit(1);}
			N_LL = std::atof(argv[2]); // Loose-Loose control region counts
			N_TL = std::atof(argv[3]); // Tight-Loose control region counts
			eff = std::atof(argv[4]); // fake rate  (must  be < 1)
			eff_err = std::atof(argv[5]); // 20% uncertainty on the measured fake rate, assuming LogNormal pdf

			cout<<"  Loose-Loose control region observation = "<<N_LL<<endl;
			cout<<"  Tight-Loose control region observation = "<<N_TL<<endl;
			cout<<"  Tight to Loose ratio = "<<eff<<endl;
			cout<<"  Uncertainty on the ratio = "<<eff_err*100<<"%"<<endl;
		}
		else{cout<<"Currently only support 1 bin with only one category of fakable object"<<endl; exit(1);}
	} 
	else {cout<<" Usage:  ./fakerate.exe 1 N_LL N_TL FakeRate FractionError"<<endl; exit(1);}

	int debug = 0;
	TMinuit *myMinuit = 0;

	int npars = 3;  // ( POI, eff, true FakeFake events	
	if(NumberOfKindsOfFakeObjects==2) npars=5;

	myMinuit = new TMinuit(npars+2);  
	if(NumberOfKindsOfFakeObjects==1)
		myMinuit->SetFCN(Function_A);
	else if (NumberOfKindsOfFakeObjects==2)
		myMinuit->SetFCN(Function_B);

	Double_t arglist[10];
	Int_t ierflg = 0;

	if(debug<100){
		arglist[0]=-1;
		myMinuit -> mnexcm("SET PRINT", arglist, 1, ierflg);
		myMinuit -> mnexcm("SET NOW", arglist, 1, ierflg);

	}

	arglist[0] = 2; // set to be 0 (was default 1), reduce lots of function calls, ---> speed up by a factor of 6 
	if(debug)cout<<" SET STRATEGY "<<arglist[0]<<endl;
	myMinuit->mnexcm("SET STRATEGY", arglist ,1,ierflg);
	//myMinuit -> mnexcm("SET NOG", arglist, 1, ierflg); // no gradiants required of FCN

	if(NumberOfKindsOfFakeObjects == 1){
		myMinuit->mnparm(0, "R0", 2., 0.01, 0, 10000, ierflg);

		//double r1 = N_TL - N_LL/(1-eff)*2*eff; if(r1<0) r1=0;
		//myMinuit->mnparm(0, "r1", r1, 0.01, 0, 100000, ierflg);
		//myMinuit->FixParameter(0);
		myMinuit->mnparm(1, "r2", N_LL/(1-eff)/(1-eff), 0.01, 0,  100000, ierflg);
		myMinuit->mnparm(2, "eff_err", 0, 0.01, -5, 5,  ierflg);
	}else if(NumberOfKindsOfFakeObjects==2){
		myMinuit->mnparm(0, "R0", 2., 0.01, 0, 10000, ierflg);
		myMinuit->mnparm(1, "r2", N_LL/(1-effe)/(1-effm), 0.01, 0,  100000, ierflg);
		myMinuit->mnparm(2, "reLmT", 10, 0.01, 0, 10000,  ierflg);
		myMinuit->mnparm(3, "effe_err", 0, 0.01, -5, 5,  ierflg);
		myMinuit->mnparm(4, "effm_err", 0, 0.01, -5, 5,  ierflg);
	}

	arglist[0] = 1;
	myMinuit->mnexcm("SET ERR", arglist ,1,ierflg);

	arglist[0] = 50000;// maximum function calls
	arglist[1] = 0.001;// tolerance
	myMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
	myMinuit->mnexcm("MINOS", arglist , 2, ierflg);


	double r, er;
	myMinuit->GetParameter(0, r, er);
	if(debug)cout<<"DELETEME before calc error,  r="<<r<<"+/-"<<er<<endl;

	// Print results
	Double_t amin,edm,errdef;
	Int_t nvpar,nparx,icstat;
	Double_t errUp, errLow, errParab=0, gcor=0; 
	myMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
	myMinuit->mnerrs(0, errUp, errLow, errParab, gcor);
	if(debug)cout<<"DELETEME errUp="<<errUp<<" errLow="<<errLow<<" errParab="<<errParab<<" gcor="<<gcor<<endl;

	double l = myMinuit->fAmin;
	myMinuit->GetParameter(0, r, er);
	if(debug)cout<<"DELETEME r="<<r<<"+/-"<<er<<" fMin="<<l<<endl;
	if(errUp==0 and errLow==0) {
		errUp = r; //er*10;
		errLow =-r; // -er*10;
	}

	if(errLow == 0 and r>0){
		errLow = -r; //-er;	
	}

	cout<<" Fitted fake contribution in Tight-Tight region:  "<<r<<" , 68\% range: [ "<<r+errLow<<"  "<<r+errUp<<" ] "<<endl;

	double R0best = r;
	double lbest = l;

	if(debug>=10) return 1;

	myMinuit->mnparm(0, "r1", r, 0.01, 0, 100000, ierflg);
	myMinuit->FixParameter(0);
	myMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
	double R0;
	double r2, er2;
	double geff, ee;
	myMinuit->GetParameter(1, r2, er2);
	myMinuit->GetParameter(2, geff, ee);
	float effcor = eff * TMath::Power(1+eff_err, geff);

	R0 = r*effcor + r2*effcor*effcor;
	if(debug)cout<<"************************* ++++++++++++++++++++"<<endl;
	if(debug)cout<< r <<" "<<r2<<" "<<effcor<<endl;


	double R2 = r2 * (1- effcor) * (1-effcor);
	double r1 =( R0best - r2 * effcor *effcor ) /effcor;
	double R1 = r2 * (1-effcor)*effcor*2 + r1 * (1-effcor); 
	double R1_ff = r2 * (1-effcor)*effcor*2;
	double R1_fp = r1 * (1-effcor); 
	double R0_ff = r2 * effcor *effcor;
	double R0_fp = r1 * effcor; // *effcor;
	double effcor_min = effcor;
	if(debug)cout << "KEEP r1="<<r<<", R0="<<R0<<" -2lnL="<< myMinuit->fAmin<<endl;

	myMinuit->mnparm(0, "r1", 0.0, 0.01, 0, 100000, ierflg);
	myMinuit->FixParameter(0);
	myMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
	myMinuit->GetParameter(1, r2, er2);
	myMinuit->GetParameter(2, geff, ee);
	effcor = eff * TMath::Power(1+eff_err, geff);

	R0 = r*effcor + r2*effcor*effcor;
	if(debug)cout << "KEEP r1="<<0<<", R0="<<R0<<" -2lnL="<< myMinuit->fAmin<<endl;


	vector<double> vR0, vl;

	double step = (errUp - errLow)/1000.; if (step <=0 ) step = 0.01;
	cout<<" step = "<<step<<endl;
	cout<<R0best<<" "<<errLow<<" "<<errUp<<endl;
	for(r=R0best+2*errLow; r<R0best+2*errUp; r+=step){
		if(r<0) continue; 
		myMinuit->mnparm(0, "r1", r, 0.01, 0, 100000, ierflg);
		myMinuit->FixParameter(0);
		myMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
		myMinuit->GetParameter(1, r2, er2);
		myMinuit->GetParameter(2, geff, ee);
		effcor = eff * TMath::Power(1+eff_err, geff);

		R0 = r*effcor + r2*effcor*effcor;
		if(debug)cout << "KEEP r1="<<r<<", R0="<<R0<<" -2lnL="<< myMinuit->fAmin<<"     effcor="<<effcor<<"   r2="<<r2<<endl;
		l = myMinuit->fAmin;
		//cout<<r<<" "<<l<<"  lbest="<<lbest<<endl;
		if(l-lbest > 2.) continue;
		vR0.push_back(r);
		vl.push_back(l);
	}
//	vR0.push_back(R0best);
//	vl.push_back(lbest);

	// plot

	TPaveText *pt = SetTPaveText(0.67, 0.4, 0.8, 0.7);
	DrawEvolution2D d2d(vR0, vl, "; fake contribution in signal region ; -2lnL ", "scanned_R0_vs_L", pt);
	d2d.draw();
	d2d.getGraph()->Draw("APC");
	double x1, x2, y1, y2;
	d2d.getCanvas()->GetRangeAxis(x1, y1, x2, y2);
	TArrow *arrowBest = new TArrow(R0best, y1+(y2-y1)/8., R0best, y1, 0.05, ">");
	arrowBest->Draw();

	double xL = interpolaton(d2d.getGraph(), lbest+1, R0best, "left");
	double xR = interpolaton(d2d.getGraph(), lbest+1, R0best, "right");
	cout<<" From plot:  fake contribution in Tight-Tight region:  "<<R0best<<" , 68\% range: [ "<<xL<<"  "<<xR<<" ] "<<endl;
	cout<<" R2 = "<<R2<<", R1="<<R1<<", effcor = "<<effcor_min<<endl;
	cout<<" R1_ff " <<R1_ff<<", R1_fp="<<R1_fp<<endl;
	cout<<" R0_ff " <<R0_ff<<", R0_fp="<<R0_fp<<endl;


	TLine *lineL = new TLine(xL, y1, xL, lbest+1);
	lineL->SetLineStyle(kDashed);
	lineL->Draw();
	TLine *lineR = new TLine(xR, y1, xR, lbest+1);
	lineR->SetLineStyle(kDashed);
	lineR->Draw();
	TLine *lineH = new TLine(xL, lbest+1, xR, lbest+1);
	lineH->SetLineStyle(kDashed);
	lineH->Draw();

	//d2d.getGraph()->Print();

	d2d.save();

}

double interpolaton(TGraph *g, double y, double xbest, TString side ){
	if(side!="left" and side!="right"){ exit(1); cout<<"fakerate::interpolaton unknown side: "<<side<<endl; }
	double x;
	g->Sort();
	int npoint = g->GetN();
	double *vx = g->GetX();
	double *vy = g->GetY();

	double xl, yl, xr, yr;
	bool byl=false, byr=false;
	for(int i=0; i<npoint; i++){
		if(side=="left") { if(vx[i]>xbest) continue;}
		else {if(vx[i]<xbest) continue;}
		if(vy[i]>=y) {
			if(byl){ if(vy[i]<yl) {yl=vy[i]; xl=vx[i];}} 
			else {yl=vy[i]; xl=vx[i]; byl=true;}
		}
		if(vy[i]<=y) {
			if(byr){ if(vy[i] > yl) {yr=vy[i]; xr=vx[i];}} 
			else {yr=vy[i]; xr=vx[i]; byr=true;}
		}
	}	
	if(!byl and !byr) return xbest;
	if(!byr) return xl; 
	if(!byl) return xr; 
	x = LinearInterpolation(xr, yr, xl, yl, y);
	return x;
}
