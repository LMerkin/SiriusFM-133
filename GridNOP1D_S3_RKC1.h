// vim:ts=2:et
#pragma  once
#include "IRP.h"
#include "Option.h"

//============================================================================//
//                          "GridNOP1D_S3_RKC1.h"                             //
// Grid Pricer for Non-IR Options, for 1D Diffusions, using 3-Point Stencils  //
//             and 1st-Order Runge-Kutta-Chebyshev Time Marshalling           //
//============================================================================//
// Non-IR means that:
// (1) Low S bound is always 0;
// (2) RateA, RateB do not depend on S;
// (3) Not suitable for very long TTE (no bounding box scaling with time).
//
namespace SiriusFM
{
  //==========================================================================//
  // "GridNOP1D_S3_RKC1" Class:                                               //
  //==========================================================================//
  template
  <
    typename Diffusion1D, typename AProvider, typename BProvider,
    typename AssetClassA, typename AssetClassB
  >
  class GridNOP1D_S3_RKC1
  {
    //------------------------------------------------------------------------//
    // Data Flds:                                                             //
    //------------------------------------------------------------------------//
  private:
    AProvider     m_irpA;
    BProvider     m_irpB;
    long          m_maxN;     // Max number of S points
    long          m_maxM;     // Max number of t points
    double* const m_grid;     // 2D Grid as a 1D array
    double* const m_S;        // S-Line
    double* const m_ts;       // TimeLine
    double* const m_ES;       // E  [S](t)
    double* const m_VarS;     // Var[S](t)  (estimated)
    int           m_N;        // Actual #of S pts
    int           m_i0;       // S[i0] = S0
    int           m_M;        // ACtual #of t pts
    bool          m_isFwd;    // Last run was Fwd?

  public:
    //------------------------------------------------------------------------//
    // Non-Default Ctor, Dtor:                                                //
    //------------------------------------------------------------------------//
    GridNOP1D_S3_RKC1
    (
      char const* a_ratesFileA,
      char const* a_ratesFileB,
      long a_maxN = 2048,
      long a_maxM = 210384
    )
    : m_irpA  (a_ratesFileA),
      m_irpB  (a_ratesFileB),
      m_maxN  (a_maxN),
      m_maxM  (a_maxM),
      m_grid  (new double[m_maxN * m_maxM]),
      m_S     (new double[m_maxN]),
      m_ts    (new double[m_maxM]),
      m_ES    (new double[m_maxM]),
      m_VarS  (new double[m_maxM]),
      m_N     (0),
      m_i0    (0),
      m_M     (0),
      m_isFwd (false)
    {
      // Zero-out all arrays:
      memset(m_grid, 0, m_maxN * m_maxM * sizeof(double));
      memset(m_S,    0, m_maxN          * sizeof(double));
      memset(m_ts,   0, m_maxM          * sizeof(double));
      memset(m_ts,   0, m_maxM          * sizeof(double));
      memset(m_VarS, 0, m_maxM          * sizeof(double));
    }

    ~GridNOP1D_S3_RKC1()
    {
      delete[] (m_grid);
      delete[] (m_S);
      delete[] (m_ts);
      delete[] (m_ES);
      delete[] (m_VarS);
      const_cast<double*&>(m_grid) = nullptr;
      const_cast<double*&>(m_S)    = nullptr;
      const_cast<double*&>(m_ts)   = nullptr;
      const_cast<double*&>(m_ES)   = nullptr;
      const_cast<double*&>(m_VarS) = nullptr;
    }

    //------------------------------------------------------------------------//
    // "Run": Performs Backward or Forward Induction:                         //
    //------------------------------------------------------------------------//
    template<bool IsFwd>
    void Run
    (
      Option<AssetClassA, AssetClassB> const*
                          a_option,         // Option Spec
      Diffusion1D const*  a_diff,
      // Grid Params:
      double              a_S0,             // S(t0); may differ from Diffusion
      time_t              a_t0,             // Abs Pricing Time
      long                a_Nints   = 500,  // #S intervals
      int                 a_tauMins = 30,   // TimeStep in minutes
      double              a_BFactor = 4.5   // #StdDevs for Upper Boundary
    );

    //------------------------------------------------------------------------//
    // "GetPxDeltaGamma0": Px, Delta and Gamma at t=0:                        //
    //------------------------------------------------------------------------//
    std::tuple<double, double, double> GetPxDeltaGamma0() const;
  };
}
