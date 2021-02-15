// vim:ts=2:et
#pragma  once
#include "GridNOP1D_S3_RKC1.h"
#include "Time.h"
#include <stdexcept>
#include <cmath>
#include <tuple>

//============================================================================//
//                         "GridNOP1D_S3_RKC1.hpp"                            //
// Grid Pricer for Non-IR Options, for 1D Diffusions, using 3-Point Stencils  //
//             and 1st-Order Runge-Kutta-Chebyshev Time Marshalling           //
//============================================================================//
namespace SiriusFM
{
  //==========================================================================//
  // "RunBI":                                                                 //
  //==========================================================================//
  template
  <
    typename Diffusion1D, typename AProvider, typename BProvider,
    typename AssetClassA, typename AssetClassB
  >
  void GridNOP1D_S3_RKC1<Diffusion1D, AProvider, BProvider,
                         AssetClassA, AssetClassB>::
  RunBI
  (
    // Option and Pricing Params:
    Option<AssetClassA, AssetClassB> const* a_option,   // Option Spec
    Diffusion1D const*  a_diff,
    // Grid Params:
    double              a_S0,        // S(t0); may differ from Diffusion
    time_t              a_t0,        // Abs Pricing Time
    long                a_Nints,     // #S intervals
    int                 a_tauMins,   // TimeStep in minutes
    double              a_BFactor    // #StdDevs for Upper Boundary
  )
  {
    //------------------------------------------------------------------------//
    // Construct the Grid:                                                    //
    //------------------------------------------------------------------------//
    assert(a_option  != nullptr && a_diff != nullptr && a_tauMins > 0 &&
           a_Nints > 0  && a_BFactor >  0);
std::cerr << "IsAmerican=" << a_option->m_isAmerican << std::endl;
    if (a_option->m_isAsian)
      throw std::invalid_argument("Asian options are not supported by 1D Grid");

    // As this is a Non-IR process, S0 must be positive:
    assert(a_S0 > 0);

    // Time to Option Expiration as Year Frac:
    double TTE = YearFracInt(a_option->m_expirTime - a_t0);

    // Fill in the TimeLine:
    // Number of t intervals:
    long tauSecs = a_tauMins * 60;
    long Mints   = (a_option->m_expirTime - a_t0) / tauSecs;

    if (TTE <= 0 || Mints <= 0)
      throw std::invalid_argument("Option has already expired, or too close");

    // Number of t points:
    m_M = Mints + 1;
    if (m_M >= m_maxM)
      throw std::invalid_argument("Too many t points");

    // Time Step:
    double tau = TTE / double(Mints);

    double integrAB = 0;
    m_ES  [0]       = a_S0;
    m_VarS[0]       = 0;

    for (int j = 0; j < m_M; ++j)
    {
      // Advance the Time Line:
      double t = YearFrac(a_t0 + j * tauSecs);
      m_ts[j]  = t;

      // Take rB(t)-rA(t) and cut the negative values, to make sure the grid
      // upper boundary is expanding with time:
      double rA = m_irpA.r(a_option->m_assetA, t);
      double rB = m_irpB.r(a_option->m_assetB, t);
      double rateDiff  = std::max<double>(rB - rA, 0);

      // Integrated rates:
      if (j < m_M-1)
      {
        integrAB    += rateDiff * tau;
        // Expected St (XXX: Do we need to store them all)?
        m_ES  [j+1]  = a_S0 * exp(integrAB);
        double sigma = a_diff->sigma(m_ES[j], t);
        m_VarS[j+1]  = m_VarS[j] + sigma * sigma * tau;
      }
    }
    // Upper Bound for S (the Lower Bound is 0):
    double B = m_ES[m_M-1] + a_BFactor * sqrt(m_VarS[m_M-1]);

    // Generate the S-Line:
    double h = B / double(a_Nints);

    // Make sure S0 is positioned exactly on the grid:
    m_i0 = int(round(a_S0 / h));  // i idx corresp to S0
    h = a_S0 / double(m_i0);
    if (!std::isfinite(h))
      throw std::invalid_argument("S0 too small, try increasing Nints");
    // Adjust the Upper Bound:
    B = h * double(a_Nints);

    m_N = a_Nints + 1;   // Number of S points
    if  (m_N > m_maxN)
      throw std::invalid_argument("Nints too large");

      // NB: The Grid is stored by-column (S-contiguous) for better locality:
    double* payOff = m_grid + (m_M-1)*m_N;

    for (int i = 0; i < m_N; ++i)
    {
      // S-line:
      m_S[i] = double(i) * h;  // Again, Low Bounds is 0

      // Create the PayOff at t=T in the Grid. Emulate "1-element" path:
      payOff[i] = a_option->Payoff(1, m_S + i, m_ts + (m_M-1));
    }
    // At Low Bound (a=0), we always have a constant boundary cond, continuous
    // with payoff:
    double fa = payOff[0];

    // At the Upper Bound, we use a const boundary cond if it is 0, otherwise
    // we fix the df/dS (a Neumann-type cond):
    bool   isNeumann = (payOff[m_N-1] != 0);
    double UBC       = isNeumann ? (payOff[m_N-1] - payOff[m_N-2]) : 0;
    for (int j = 0; j < m_M-1; ++j)
      m_grid[j*m_N] = fa;  // Low Bound

    // Time Marshalling:
    double D2 = 2 * h * h;     // Denom in the Diffusive Term

    for (int j = m_M-1; j >= 1; --j)
    {
      double const* fj   = m_grid + j * m_N;  // Prev time layer (j)
      double*       fjm1 = const_cast<double*>(fj - m_N);
                    // Curr time layer to be filled in (j-1)

      double tj      = m_ts[j];
      double rateAj  = m_irpA.r(a_option->m_assetA, tj);
      double rateBj  = m_irpB.r(a_option->m_assetB, tj);
      double C1      = (rateBj - rateAj) / (2 * h); // Coeff in the Conv Term

      fjm1[0]   = fa;   // Low Bound

//#     pragma acc parallel loop copyin(fj[0:m_N],fjm1[0:m_N-1],tj,C1) copyout(fjm1[0:m_N-1])
#     pragma omp parallel for
      for (int i = 1; i <= m_N-2; ++i)
      {
        double Si    = m_S [i];
        double fjim1 = fj[i-1];
        double fji   = fj[i];
        double fjip1 = fj[i+1];
        double sigma = a_diff->sigma(Si, tj);

        double DfDt  = rateBj  * fji              // Reactive Term
                     - C1 * Si * (fjip1 - fjim1)  // Conv Term
                     - sigma * sigma / D2 * (fjip1 - 2*fji + fjim1);
        // FIXME: Euler's method instead of RKC1:
        fjm1[i]      = fji - tau * DfDt;
      }
      fjm1[m_N-1] = isNeumann ? (fjm1[m_N-2] + UBC) : UBC;

      // Grid allows us to price American options as well:
      if (a_option->m_isAmerican)
        for (int i = 0; i < m_N; ++i)
        {
          // Intrinsic value of the option is the payoff evaluated under the
          // curr underlying px Si:
          double  intrVal = a_option->Payoff(1, m_S + i, &tj);
/*
if (intrVal > fjm1[i])
  std::cerr << "t=" << tj << ", S=" << m_S[i] << ": IntrVal=" << intrVal << ", OptPx=" << fjm1[i] << std::endl;
*/
          fjm1[i] = std::max<double>(fjm1[i],  intrVal);
        }
    }
    // End of Time Marshalling
  }

  //==========================================================================//
  // "GetPriceDeltaGamma0": Assume S0 is on the Grid:                         //
  //==========================================================================//
  template
  <
    typename Diffusion1D, typename AProvider, typename BProvider,
    typename AssetClassA, typename AssetClassB
  >
  std::tuple<double, double, double>
  GridNOP1D_S3_RKC1<Diffusion1D, AProvider, BProvider,
                    AssetClassA, AssetClassB>::GetPxDeltaGamma0() const
  {
    if (m_M == 0 || m_N == 0)
      throw std::runtime_error("RunBI first!");

    assert(0 <= m_i0 && m_i0 < m_N);

    double h     = m_S[1] - m_S[0];
    double px    = m_grid[m_i0];    // j=0
    double delta = 0;
    double gamma = 0;
    if (0 < m_i0 && m_i0 <= m_N-2)
    {
      delta = (m_grid[m_i0+1] - m_grid[m_i0-1]) / (2*h);
      gamma = (m_grid[m_i0+1] - 2*m_grid[m_i0] + m_grid[m_i0-1]) / (h*h);
    }
    else
    if (m_i0 == 0)
      delta = (m_grid[1]   - m_grid[0])   / h;    // gamma remains 0
    else
    {
      assert(m_i0  == m_N-1);
      delta = (m_grid[m_N-1] - m_grid[m_N-2]) / h; // gamma remains 0
    }
    return std::make_tuple(px, delta, gamma);
  }
}
