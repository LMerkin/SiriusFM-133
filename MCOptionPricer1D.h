#pragma  once
#include "MCEngine1D.h"
#include "Option.h"

namespace SiriusFM
{
  //=========================================================================//
  // "MCOptionPricer1D":                                                     //
  //=========================================================================//
  template
  <
    typename Diffusion1D, typename AProvider, typename BProvider,
    typename AssetClassA, typename AssetClassB
  >
  class MCOptionPricer1D
  {
  private:
    //=======================================================================//
    // Path Evaluator for Option Pricing:                                    //
    //=======================================================================//
    class OPPathEval
    {
    private:
      Option<AssetClassA, AssetClassB> const* const m_option;
      long   m_P;     // Total paths evaluated
      double m_sum;   // Sum of Payoffs
      double m_sum2;  // Sum of Payoff^2
      double m_minPO; // Min PayOff
      double m_maxPO; // Max PayOff

    public:
      OPPathEval(Option<AssetClassA, AssetClassB> const* a_option)
      : m_option(a_option),
        m_P     (0),
        m_sum   (0),
        m_sum2  (0),
        m_minPO ( INFINITY),
        m_maxPO (-INFINITY)
      
      { assert(m_option != nullptr); }

      void operator() (long a_L,     long a_PM,
                       double const* a_paths, double const* a_ts)
      {
        for (long p = 0; p < a_PM; ++p)
        {
          double const* path = a_paths + p * a_L;
          double payOff      = m_option->Payoff(a_L, path, a_ts);
          m_sum  += payOff;
          m_sum2 += payOff * payOff;
          m_minPO = std::min<double>(m_minPO, payOff);
          m_maxPO = std::max<double>(m_maxPO, payOff);
        }
        m_P += a_PM;
      }

      // GetPx: Returns E[Px]:
      //
      double GetPx() const
      {
        if (m_P < 2)
          throw std::runtime_error("Empty OPPathEval");
        return m_sum / double(m_P);
      }

      // GetStats: Returns StD[Px], Min[PayOff], Max[PayOff]:
      // (for debugging only):
      //
      std::tuple<double, double, double> GetStats() const
      {
        if (m_P < 2)
          throw std::runtime_error("Empty OPPathEval");
        double px  =  m_sum  / double(m_P);
        double var = (m_sum2 - double(m_P) * px * px) / double(m_P - 1);
        assert(var >= 0);
        return std::make_tuple(sqrt(var), m_minPO, m_maxPO);
      }
    };
    //=======================================================================//
    // Flds:                                                                 //
    //=======================================================================//
    Diffusion1D const* const  m_diff;
    AProvider                 m_irpA;
    BProvider                 m_irpB;
    MCEngine1D<Diffusion1D, AProvider, BProvider, AssetClassA, AssetClassB,
               OPPathEval>    m_mce;
    bool                      m_useTimerSeed;

  public:
    //-----------------------------------------------------------------------//
    // Non-Default Ctor:                                                     //
    //-----------------------------------------------------------------------//
    MCOptionPricer1D
    (
      Diffusion1D const*  a_diff,
      const char*         a_irsFileA,
      const char*         a_irsFileB,
      bool                a_useTimerSeed
    )
    : m_diff(a_diff),
      m_irpA(a_irsFileA),
      m_irpB(a_irsFileB),
      m_mce (102271, 4096), // (5-min points in 1y) * 4k paths in-memory
      m_useTimerSeed(a_useTimerSeed)
    {}

    //-----------------------------------------------------------------------//
    // The Pricing Function:                                                 //
    //-----------------------------------------------------------------------//
    double Px
    (
      // Instrument Spec:
      Option<AssetClassA, AssetClassB> const* a_option,
      // Pricing Time:
      time_t              a_t0,
      // MC Params:
      int                 a_tauMins = 15,
      long                a_P       = 100'000
    );
  };
}
