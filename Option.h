#pragma once 
#include <ctime>

namespace SiriusFM
{
	class Option
	{
  public:
    time_t const m_expirTime;
		bool   const m_isAmerican;

		Option(time_t  a_expirTime, bool a_isAmerican)
    : m_expirTime (a_expirTime),
      m_isAmerican(a_isAmerican)
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
}
