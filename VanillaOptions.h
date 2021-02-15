#pragma once
#include "Option.h"
#include <cmath>
#include <cassert>

namespace SiriusFM
{
  //=========================================================================//
  // Generic European or American (but not Asian) Call:                      //
  //=========================================================================//
  template<typename AssetClassA,  typename AssetClassB>
	class CallOption final: public Option<AssetClassA, AssetClassB>
	{
		double const m_K;
	public:
		CallOption
    (
      AssetClassA a_assetA,
      AssetClassB a_assetB,
      double      a_K,
      time_t      a_expirTime,
      bool        a_isAmerican
    )
    : Option<AssetClassA, AssetClassB>
            (a_assetA, a_assetB, a_expirTime, a_isAmerican, false),
            // IsAsian=false
	    m_K(a_K)
		{
			if(a_K <= 0)
				throw std::invalid_argument("Bad K");
		}

		~CallOption() override {}

		virtual double Payoff(long a_L, 
							  double const*  a_path,
							  double const*  a_ts = nullptr) const override
		{
			assert(a_L > 0 && a_path != nullptr);
			return std::max<double>(a_path[a_L - 1] - m_K, 0);
		}
	};

  //=========================================================================//
  // Generic European or American (but not Asian) Put:                       //
  //=========================================================================//
  template<typename  AssetClassA,  typename AssetClassB>
	class PutOption final: public Option<AssetClassA, AssetClassB>
	{
		double const m_K;

	public:
		PutOption
    (
      AssetClassA a_assetA,
      AssetClassB a_assetB,
      double      a_K,
      time_t      a_expirTime,
      bool        a_isAmerican
    )
    : Option<AssetClassA, AssetClassB>
            (a_assetA, a_assetB, a_expirTime, a_isAmerican, false),
            // IsAsian=false
	    m_K(a_K)
		{
			if(a_K <= 0)
				throw std::invalid_argument("Bad K");
		}

		~PutOption() override {}

		virtual double Payoff(long a_L, 
							  double const*  a_path,
							  double const*  a_ts = nullptr) const override
		{
			assert(a_L > 0 && a_path != nullptr);
			return std::max<double>(m_K - a_path[a_L - 1], 0);
		}
	};

  //=========================================================================//
  // Aliases:                                                                //
  //=========================================================================//
  using CallOptionFX = CallOption<CcyE, CcyE>;
  using PutOptionFX  = PutOption <CcyE, CcyE>;
}
