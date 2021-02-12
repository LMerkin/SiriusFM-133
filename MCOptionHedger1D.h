#pragma  once
#include "MCEngine1D.h"
#include "Option.h"
#include <functional>

namespace SiriusFM
{
  //=========================================================================//
  // "MCOptionHedger1D":                                                     //
  //=========================================================================//
  template
  <
    typename Diffusion1D, typename AProvider, typename BProvider,
    typename AssetClassA, typename AssetClassB
  >
  class MCOptionHedger1D
  {
  public:
      using DeltaFunc =
        std::function<double(double, double)>; // (S,t) -> delta

  private:
    //=======================================================================//
    // Path Evaluator for Option Hedging:                                    //
    //=======================================================================//
    class OHPathEval
    {
    private:
      Option<AssetClassA, AssetClassB> const* const m_option;
      AProvider                        const* const m_irpA;
      BProvider                        const* const m_irpB;
      double*                                       m_ratesA;
      double*                                       m_ratesB;
      // Hedging policy:
      double const           m_C0;       // Initial option premium
      DeltaFunc const* const m_DeltaFunc;
      double const           m_DeltaAcc; // Delta rounded to a multiple of this
      // Monte-Carlo stats:
      long                   m_P;        // Total paths evaluated
      double                 m_sumPnL;   // Sum of Residual P&Ls
      double                 m_sumPnL2;  // Sum of PnL^2s
      double                 m_minPnL;   // Min PnL
      double                 m_maxPnL;   // Max PnL

    public:
      OHPathEval
      (
        Option<AssetClassA, AssetClassB> const* a_option,
        AProvider const*                        a_irpA,
        BProvider const*                        a_irpB,
        double                                  a_C0,
        DeltaFunc const*                        a_deltaFunc,
        double                                  a_deltaAcc
      )
      : m_option    (a_option),
        m_irpA      (a_irpA),
        m_irpB      (a_irpB),
        m_ratesA    (nullptr),
        m_ratesB    (nullptr),
        m_C0        (a_C0),
        m_DeltaFunc (a_deltaFunc),
        m_DeltaAcc  (a_deltaAcc),
        m_P         (0),
        m_sumPnL    (0),
        m_sumPnL2   (0),
        m_minPnL    ( INFINITY),
        m_maxPnL    (-INFINITY)
      {
        assert(m_option != nullptr && m_DeltaFunc != nullptr &&
               m_DeltaAcc >= 0.0   && m_irpA      != nullptr &&
               m_irpB   != nullptr);
      }

      // Dtor:
      ~OHPathEval()
      {
        delete[] (m_ratesA);
        delete[] (m_ratesB);
        m_ratesA = nullptr;
        m_ratesB = nullptr;
      }

      void operator() (long a_L,     long a_PM,
                       double const* a_paths, double const* a_ts)
      {
        // If rates are not yet available, pre-compute them:
        if (m_ratesA == nullptr)
        {
          m_ratesA = new double[a_L];
          for (long l = 0; l < a_L; ++l)
            m_ratesA[l] = m_irpA->r(m_option->m_assetA, a_ts[l]);
        }

        if (m_ratesB == nullptr)
        {
          m_ratesB = new double[a_L];
          for (long l = 0; l < a_L; ++l)
            m_ratesB[l] = m_irpB->r(m_option->m_assetB, a_ts[l]);
        }

        // Evaluate all stored paths:
        for (long p = 0; p < a_PM; ++p)
        {
          double const* path = a_paths + p * a_L;

          // Perform Delta-hedging along this path:
          double M     = - m_C0;   // We long the option, short C0: curr Money
          double delta = 0.0;      // Curr delta

          for (long l = 0; l < a_L; ++l)
          {
            double St     = path[l];   // Curr underlying px
            double t      = a_ts[l];   // Curr time

            if (l > 0)
            {
              // Manage the money account:
              double tau = t - a_ts[l-1];
              double Sp  = path[l-1];
              M += M  * tau * m_ratesB[l-1];

              // Also dividends (wrt prev S):
              M += Sp * tau * m_ratesA[l-1];
            }

            // Delta Hedging (no need for it in the last point):
            if (l < a_L - 1)
            {
              double deltaN = (*m_DeltaFunc)(St, t);
              // Round "deltaN" to a multiple of "DeltaAcc":
              // Also,  deltaN changes sign (as we long the option)
              deltaN = - round(deltaN / m_DeltaAcc) * m_DeltaAcc;

              if (delta != deltaN)
              {
                // Re-Hedge:
                M -= (deltaN - delta) * St;
                delta = deltaN;
              }
            }
          }
          // End of Path. Get the PayOff and the total portfolio value:
          double PnL =
            M + delta * path[a_L-1] + m_option->Payoff(a_L, path, a_ts);
// std::cerr << "M=" << M << ", delta=" << delta    << ", payoff="
//        << (m_option->Payoff(a_L, path, a_ts)) << std::endl;

          // Update the stats:
          m_sumPnL  += PnL;
          m_sumPnL2 += PnL * PnL;
          m_minPnL  = std::min<double>(m_minPnL, PnL);
          m_maxPnL  = std::max<double>(m_maxPnL, PnL);
        }
        // Increment the number of paths processed:
        m_P += a_PM;
      }

      // GetStats: Returns E[PnL], StD[PnL], Min[PnL], Max[PnL]:
      //
      std::tuple<double, double, double, double> GetStats() const
      {
        if (m_P < 2)
          throw std::runtime_error("Empty OHPathEval");
        double mean =  m_sumPnL  / double(m_P);
        double var  = (m_sumPnL2 - double(m_P) * mean * mean)
                      / double(m_P - 1);
        assert(var >= 0);
        return std::make_tuple(mean, sqrt(var), m_minPnL, m_maxPnL);
      }
    };
    //=======================================================================//
    // Flds:                                                                 //
    //=======================================================================//
    Diffusion1D const* const  m_diff;
    AProvider                 m_irpA;
    BProvider                 m_irpB;
    MCEngine1D<Diffusion1D, AProvider, BProvider, AssetClassA, AssetClassB,
               OHPathEval>    m_mce;
    bool                      m_useTimerSeed;

  public:
    //-----------------------------------------------------------------------//
    // Non-Default Ctor:                                                     //
    //-----------------------------------------------------------------------//
    MCOptionHedger1D
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
    // The Hedging Simulator:                                                //
    //-----------------------------------------------------------------------//
    std::tuple<double, double, double, double> SimulateHedging
    (
      Option<AssetClassA, AssetClassB> const* a_option,
      time_t              a_t0,               // Start time
      double              a_C0,
      DeltaFunc const*    a_DeltaFunc,
      double              a_DeltaAcc,
      int                 a_tauMins = 15,
      long                a_P       = 100'000
    );

    //-----------------------------------------------------------------------//
    // Accessors:                                                            //
    //-----------------------------------------------------------------------//
    double GetRateA(AssetClassA a_assetA, double a_ty) const
      { return m_irpA.r(a_assetA, a_ty); }

    double GetRateB(AssetClassB a_assetB, double a_ty) const
      { return m_irpB.r(a_assetB, a_ty); }
  };
}
