#ifndef  BAYESIANLIMIT_H
#define  BAYESIANLIMIT_H
#include "CRandom.h"
#include <string>
#include <vector>
#include "CountingModel.h"
using namespace std;
namespace lands{
class BayesianLimit
{
	
	public:
	BayesianLimit(CountingModel *cms, int nexps, CRandom *rdm);
	~BayesianLimit();
	void SetRdm(CRandom *rnd);
	void SetDebug(int d){debug=d;};
	void Run();
	double Limit_mean();
	double Limit_sigma(int nsigma);
	double LimitForExpectedBonly(); //Asimov 

	vector<double> GetDifferentialLimits();
	vector<double> GetvLimits();
	vector<double> GetvLimits_CP();

	private:
	//--- input
	int _nexps;
	CRandom *_rdmGen;
	int _nchannels;
	double _nsig[10000]; // "maximum 10000 channels", FIXME
	double _nbkg[10000];

	//--- return values	
	double _meanLimit;
	double _medianLimit;
	double _1SigmaLow;
	double _1SigmaHigh;
	double _2SigmaLow;
	double _2SigmaHigh;

	vector<double> v_rLimits;
	vector<int> v_rLimits_n;
	vector<double> v_rLimits_sorted;
	vector<double> v_cumulaP_sorted;
	vector<double> _differentialR95s;

	void SingleChannelRun(double nsig, double nbkg);

	int debug;
	CountingModel *_cms;
};
};
#endif   /* ----- #ifndef BAYESIANLIMIT_H----- */

