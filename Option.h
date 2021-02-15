#pragma once 
#include "IRP.h"
#include <ctime>

namespace SiriusFM
{
  //=========================================================================//
  // Fully-Generic "Option":                                                 //
  //=========================================================================//
  template<typename AssetClassA, typename AssetClassB>
	class Option
	{
  public:
    AssetClassA const m_assetA;  // Option Underlying Instrument: A/B
    AssetClassB const m_assetB;
    time_t      const m_expirTime;
		bool        const m_isAmerican;
    bool        const m_isAsian;

		Option
    (
      AssetClassA a_assetA,
      AssetClassB a_assetB,
      time_t      a_expirTime,
      bool        a_isAmerican,
      bool        a_isAsian
    )
    : m_assetA    (a_assetA),
      m_assetB    (a_assetB),
      m_expirTime (a_expirTime),
      m_isAmerican(a_isAmerican),
      m_isAsian   (a_isAsian)
    {}

		virtual double Payoff
    (
      long a_L, 
		  double const* a_path,
			double const* a_ts
    )
    const = 0;

		virtual ~Option() {};
	};

  //=========================================================================//
  // Alias: "OptionFX":                                                      //
  //=========================================================================//
  using OptionFX = Option<CcyE, CcyE>;
}
