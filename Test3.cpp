#include "IRProviderConst.h"
#include "DiffusionGBM.h"
#include "VanillaOptions.h"
#include "MCOptionPricer1D.hpp"
#include <iostream>
#include <cstring>

using namespace SiriusFM;
using namespace std;

int main(int argc, char** argv)
{
	if(argc != 9)
	{
		cerr << "params: mu, sigma, S0,\nCall/Put, K, Tdays,\ntau_mins, P\n";
		return 1;
	}
	double mu = atof(argv[1]);
	double sigma = atof(argv[2]);
	double S0 = atof(argv[3]);
	const char* OptType = argv[4];
	double K = atof(argv[5]);
	long T_days = atol(argv[6]);
	int tau_mins = atoi(argv[7]);
	long P = atol(argv[8]);

	assert(sigma > 0 &&
		   S0 > 0 &&
		   T_days > 0 &&
		   tau_mins > 0 &&
		   P > 0 &&
		   K > 0);

	CcyE ccyA = CcyE::USD;
	CcyE ccyB = CcyE::USD;

  char const* ratesFileA = nullptr;
  char const* ratesFileB = nullptr;
  bool useTimerSeed = true;

	DiffusionGBM diff(mu, sigma, S0);

  // The following Pricer is for FX (CcyE / CcyE):
  MCOptionPricer1D<decltype(diff), IRPConst, IRPConst, CcyE, CcyE>
    pricer(&diff, ratesFileA, ratesFileB, useTimerSeed);

  // Create the Option spec:
	time_t t0 = time(nullptr);   // Pricing Time
	time_t T  = t0 + SEC_IN_DAY * T_days;

	OptionFX const* opt = nullptr;
  if (strcmp(OptType, "Call") == 0)
		opt = new EurCallOptionFX(ccyA, ccyB, K, T);
  else
	if (strcmp(OptType, "Put") == 0)
	  opt = new EurPutOptionFX (ccyA, ccyB, K, T);
  else
		throw invalid_argument("Bad option type");

  // Presto! Run the Pricer:
  double px = pricer.Px(opt, t0, tau_mins, P);

  cout << "Px=" << px << endl;
  delete opt;
	return 0;
}
