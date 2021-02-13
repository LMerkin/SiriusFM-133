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
    : m_irpA(a_ratesFileA),
      m_irpB(a_ratesFileB),
      m_maxN(a_maxN),
      m_maxM(a_maxM),
      m_grid(new double[m_maxN * m_maxM]),
      m_S   (new double[m_maxN]),
      m_ts  (new double[m_maxM]),
      m_ES  (new double[m_maxM]),
      m_VarS(new double[m_maxM])
    {}

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
    // "RunBI": Performs Backward-Induction:                                  //
    //------------------------------------------------------------------------//
    void RunBI
    (
      Option<AssetClassA, AssetClassB> const* a_option,   // Option Spec
      Diffusion1D const*  a_diff,
      // Grid Params:
      double              a_S0,             // S(t0); may differ from Diffusion
      time_t              a_t0,             // Abs Pricing Time
      long                a_N       = 500,  // #S points
      int                 a_tauMins = 30,   // TimeStep in minutes
      double              a_BFactor = 4.5   // #StdDevs for Upper Boundary
    );
  };
}
