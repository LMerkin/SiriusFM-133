#include "DiffusionGBM.h"
#include "IRProviderConst.h"
#include "MCEngine1D.hpp"

using namespace SiriusFM;
using namespace std;

int main(int argc, char** argv)
{
	if(argc != 7)
	{
		cerr << "params: mu, sigma, S0, T_days, tau_mins, P\n";
		return 1;
	}
	double mu = atof(argv[1]);
	double sigma = atof(argv[2]);
	double S0 = atof(argv[3]);
	long T_days = atol(argv[4]);
	int tau_mins = atoi(argv[5]);
	long P = atol(argv[6]);

	//check sigma > 0, S0 > 0, T_days > 0, tau_min > 0, P > 0

	CcyE ccyA = CcyE::USD;
	CcyE ccyB = CcyE::USD;

	IRProvider<IRModeE::Const> irp(nullptr);
	DiffusionGBM diff(mu, sigma);

	MCEngine1D<DiffusionGBM, decltype(irp),
		decltype(irp), CcyE, CcyE> mce(20000, 20000);

	time_t t0 = time(nullptr);
	cerr << "t0 = " << t0 << endl;
	time_t T = t0 + SEC_IN_DAY * T_days;
	cerr << "T = " << T << endl;
	double Ty = double(T_days)/AVG_DAYS_IN_YEAR;
	cerr << "Ty = " << Ty << endl;

	//Run MC
	cerr << "Running Simulate() with: " << endl;
	cerr << t0 << endl << T << endl << tau_mins << endl << P << endl << S0;
	cerr << endl;
	mce.Simulate<false>(t0, T, tau_mins, P, S0, &diff, &irp, &irp, ccyA, ccyA);
	//S0 better be given to diff, not to MCE

	//Analyse the result
	auto res = mce.GetPaths();
	long L1 = get<0>(res);
	long P1 = get<1>(res);
	double const* paths = get<2>(res);

	//compute E of log S_T
	double EST = 0.0;
	double EST2 = 0.0;
	int N = 0; //valid paths
	for(int p = 0; p < P1; ++p)
	{
		double const* path = paths + p * L1;
		double ST = path[L1-1];
		//in pratice may get ST <= 0
		if(ST <= 0)
			continue;
		++N;
		double RT = log(ST/S0);	
		EST += RT;
		EST2 += RT * RT;
	}

	assert(N > 1);
	EST /= double(N); //(mu - sigma^2/2) * T
	double VarST = (EST2 - double(N) * EST * EST )/ double(N-1); //sigma^2 * T

	double sigma2E = VarST / Ty;
	double muE = (EST + VarST / 2) / Ty;

	cout << "mu = " << mu << ", mu_est = " << muE << endl;
	cout << "sigma^2 = " << sigma * sigma << ", sigma^2_est = " << sigma2E;
	cout << endl;

	return 0;
}
