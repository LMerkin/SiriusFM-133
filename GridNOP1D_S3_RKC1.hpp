// vim:ts=2:et
#pragma  once
#include "GridNOP1D_S3_RKC1.h"
#include "Time.h"
#include <cmath>
#include <stdexcept>

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
    long                a_N,         // #points
    int                 a_tauMins,   // TimeStep in minutes
    double              a_BFactor    // #StdDevs for Upper Boundary
  )
  {
    //------------------------------------------------------------------------//
    // Construct the Grid:                                                    //
    //------------------------------------------------------------------------//
    assert(a_option != nullptr && a_diff != nullptr && a_tauMins > 0 &&
           a_N > 1  && a_BFactor >  0);
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
    long M     = Mints + 1;
    if (M >= m_maxM)
      throw std::invalid_argument("Too many t points");

    // Time Step:
    double tau = TTE / double(Mints);

    double integrAB = 0;
    m_ES  [0]       = a_S0;
    m_VarS[0]       = 0;

    for (int j = 0; j < M; ++j)
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
      if (j < M-1)
      {
        integrAB    += rateDiff * tau;
        // Expected St (XXX: Do we need to store them all)?
        m_ES  [j+1]  = a_S0 * exp(integrAB);
        double sigma = a_diff->sigma(m_ES[j], t);
        m_VarS[j+1]  = m_VarS[j] + sigma * sigma * tau;
      }
    }
    // Upper Bound for S (the Lower Bound is 0):
    double B = m_ES[M-1] + a_BFactor * sqrt(m_VarS[M-1]);

    // Generate the S-Line:
    double h = B / double(a_N-1);
  }
}
