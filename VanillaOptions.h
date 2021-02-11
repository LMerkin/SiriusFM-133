#pragma once
#include "Option.h"
#include <cmath>
#include <cassert>

namespace SiriusFM
{
  //=========================================================================//
  // Generic European Call:                                                  //
  //=========================================================================//
  template<typename AssetClassA,  typename AssetClassB>
	class EurCallOption final: public Option<AssetClassA, AssetClassB>
	{
		double const m_K;
	public:
		EurCallOption
    (
      AssetClassA a_assetA,
      AssetClassB a_assetB,
      double      a_K,
      time_t      a_expirTime
    )
    : Option<AssetClassA, AssetClassB>(a_assetA, a_assetB, a_expirTime, false),
                                                // IsAmerican=false
	    m_K(a_K)
		{
			if(a_K <= 0)
				throw std::invalid_argument("Bad K");
		}

		~EurCallOption() override {}

		virtual double Payoff(long a_L, 
							  double const*  a_path,
							  double const*  a_ts = nullptr) const override
		{
			assert(a_L > 0 && a_path != nullptr);
			return std::max<double>(a_path[a_L - 1] - m_K, 0);
		}
	};

  //=========================================================================//
  // Generic European Put:                                                   //
  //=========================================================================//
  template<typename  AssetClassA,  typename AssetClassB>
	class EurPutOption final: public Option<AssetClassA, AssetClassB>
	{
		double const m_K;

	public:
		EurPutOption
    (
      AssetClassA a_assetA,
      AssetClassB a_assetB,
      double      a_K,
      time_t      a_expirTime
    )
    : Option<AssetClassA, AssetClassB>(a_assetA, a_assetB, a_expirTime, false),
                                                          // IsAmerican=false
	    m_K(a_K)
		{
			if(a_K <= 0)
				throw std::invalid_argument("Bad K");
		}

		~EurPutOption() override {}

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
  using EurCallOptionFX = EurCallOption<CcyE, CcyE>;
  using EurPutOptionFX  = EurPutOption <CcyE, CcyE>;
}
