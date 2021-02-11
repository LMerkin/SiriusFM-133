#pragma once
#include "MCOptionPricer1D.h"
#include "MCEngine1D.hpp"

namespace SiriusFM
{
  //=========================================================================//
  // "MCOptionPricer1D::Px":                                                 //
  //=========================================================================//
  template
  <
    typename Diffusion1D, typename AProvider, typename BProvider,
    typename AssetClassA, typename AssetClassB
  >
  double MCOptionPricer1D<Diffusion1D, AProvider, BProvider,
                          AssetClassA, AssetClassB>::
  Px
  (
    // Instrument Spec:
    Option<AssetClassA, AssetClassB> const* a_option,
    // Pricing Time:
    time_t              a_t0,
    // MC Params:
    int                 a_tauMins,
    long                a_P
  )
  {
    assert(a_option != nullptr && a_tauMins > 0 && a_P > 0);

    // Path Evaluator:
    OPPathEval pathEval(a_option);

	  // Run MC: Option pricing is Risk-Neutral:
	  m_mce.template Simulate<true>
      (a_t0,   a_option->m_expirTime, a_tauMins, a_P, m_useTimerSeed,
       m_diff, &m_irpA, &m_irpB,  a_option->m_assetA, a_option->m_assetB,
       &pathEval);

    // Get the price from PathEval:
    double px = pathEval.GetPx();

    // Apply the Discount Factor on B:
    px *= m_irpB.DF(a_option->m_assetB, a_t0, a_option->m_expirTime);
    return px;
  }
}
