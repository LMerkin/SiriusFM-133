#include "IRProviderConst.h"
#include "DiffusionGBM.h"
#include "VanillaOptions.h"
#include "MCOptionHedger1D.hpp"
#include <iostream>
#include <cstring>
#include <cmath>

using namespace SiriusFM;
using namespace std;

namespace
{
  // CDF of the Standard Normal:
  inline double Phi(double x)
    { return 0.5 * (1.0 + erf(x / M_SQRT2)); }

  // BSM Pricer (const params):
  // TTE is Time-To-Expiration (as Year Fraction):
  //
  inline double BSMPxCall(double a_S0,    double a_K,     double a_TTE,
                          double a_rateA, double a_rateB, double a_sigma)
  {
    assert(a_S0 > 0 && a_K > 0 && a_sigma > 0);
    if (a_TTE <= 0)
      // Return PayOff:
      return std::max<double>(a_S0 - a_K, 0);

    double xd = a_sigma * sqrt(a_TTE);
    double x1 =
      (log(a_S0 / a_K) +
        (a_rateB - a_rateA + a_sigma * a_sigma / 2.0) * a_TTE) / xd;
    double x2 = x1 - xd;
    double px = a_S0 * exp(-a_rateA * a_TTE) * Phi(x1) -
                a_K  * exp(-a_rateB * a_TTE) * Phi(x2);
    return px;
  }

  inline double BSMPxPut (double a_S0,    double a_K,     double a_TTE,
                          double a_rateA, double a_rateB, double a_sigma)
  {
    double px =
      BSMPxCall(a_S0, a_K, a_TTE, a_rateA, a_rateB, a_sigma) - a_S0 +
      exp(-a_rateB * a_TTE) * a_K;
    assert(px > 0.0);
    return px;
  }

  // Deltas of Call and Put:
  inline double BSMDeltaCall(double a_S0,    double a_K,     double a_TTE,
                             double a_rateA, double a_rateB, double a_sigma)
  {
    assert(a_S0 > 0 && a_K > 0 && a_sigma > 0);
    if (a_TTE <= 0)
      return (a_S0 < a_K) ? 0 : (a_S0 > a_K) ? 1 : 0.5;

    double xd = a_sigma * sqrt(a_TTE);
    double x1 =
      (log(a_S0 / a_K) +
        (a_rateB - a_rateA + a_sigma * a_sigma / 2.0) * a_TTE) / xd;
    return Phi(x1);
  }

  inline double BSMDeltaPut (double a_S0,    double a_K,     double a_TTE,
                             double a_rateA, double a_rateB, double a_sigma)
  {
    return BSMDeltaCall(a_S0, a_K, a_TTE, a_rateA, a_rateB, a_sigma) - 1.0;
  }
}

int main(int argc, char** argv)
{
	if(argc != 10)
	{
		cerr << "params: mu, sigma, S0,\nCall/Put, K, Tdays,\ndeltaAcc,\ntau_mins, P\n";
		return 1;
	}
	double mu           = atof(argv[1]);
	double sigma        = atof(argv[2]);
	double S0           = atof(argv[3]);
	const char* OptType = argv[4];
	double K            = atof(argv[5]);
	long T_days         = atol(argv[6]);
  double deltaAcc     = atof(argv[7]);
	int tau_mins        = atoi(argv[8]);
	long P              = atol(argv[9]);

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

  // The following Hedger is for FX (CcyE / CcyE):
  MCOptionHedger1D<decltype(diff), IRPConst, IRPConst, CcyE, CcyE>
    hedger(&diff, ratesFileA, ratesFileB, useTimerSeed);

  // Create the Option spec:
	time_t t0  = time(nullptr);            // Abs Start Time
	time_t T   = t0 + SEC_IN_DAY * T_days; // Abs Expir Time in Sec from Epoch
  double TTE = YearFracInt(T - t0);
  double Ty  = EPOCH_BEGIN + double(T) / AVG_SEC_IN_YEAR;
                                         // Expir Time as YYYY.YearFrac
	OptionFX const* opt = nullptr;
  decltype(hedger)::DeltaFunc const* deltaFunc = nullptr;
  double          C0  = 0.0;

  // Rates are constant here:
  double rateA = hedger.GetRateA(ccyA, 0.0);  // Any t
  double rateB = hedger.GetRateA(ccyB, 0.0);  // Any t

  // Closures for Deltas:
  function<double(double,double)> deltaCall
    (
      [K, Ty, rateA, rateB, sigma]
      (double a_St, double a_t) -> double
      {
        double currTTE = Ty - a_t;
        return BSMDeltaCall(a_St, K, currTTE, rateA, rateB, sigma);
      }
    );
  function<double(double,double)> deltaPut
    (
      [K, Ty, rateA, rateB, sigma]
      (double a_St, double a_t) -> double
      {
        double currTTE = Ty - a_t;
        return BSMDeltaPut (a_St, K, currTTE, rateA, rateB, sigma);
      }
    );

  if (strcmp(OptType, "Call") == 0)
  {
		opt       = new EurCallOptionFX(ccyA, ccyB, K, T);
    C0        = BSMPxCall(S0, K, TTE, rateA, rateB, sigma);
    deltaFunc = &deltaCall;
  }
  else
	if (strcmp(OptType, "Put") == 0)
  {
	  opt       = new EurPutOptionFX (ccyA, ccyB, K, T);
    C0        = BSMPxPut(S0, K, TTE, rateA, rateB, sigma);
    deltaFunc = &deltaPut;
  }
  else
		throw invalid_argument("Bad option type");

  // Presto! Run the Hedger:
  auto res = hedger.SimulateHedging
    (opt, t0, C0, deltaFunc, deltaAcc, tau_mins, P);

  double EPnL   = get<0>(res);
  double StDPnL = get<1>(res);
  double MinPnL = get<2>(res);
  double MaxPnL = get<3>(res);

  cout << "E[PnL] = "    << EPnL   << ", StD[PnL] = " << StDPnL
       << ", Min[PnL] =" << MinPnL << ", Max[PnL] = " << MaxPnL << endl;

  delete opt;
	return 0;
}
