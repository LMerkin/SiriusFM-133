#pragma once

#include<stdexcept>
namespace SiriusFM
{
                
	class DiffusionGBM
	{
		double const m_muBar;
		double const m_sigmaBar;
		double const m_S0;

	public:
		DiffusionGBM(double a_m, double a_s, double a_s0): m_muBar(a_m), 
													       m_sigmaBar(a_s),
													       m_S0(a_s0)
		{
			if(m_sigmaBar <= 0)
			{
				throw std::invalid_argument("Bad sigma");
			}
		};

		double mu(double a_St, double a_t) const {return m_muBar * a_St;};
		double sigma(double a_St, double a_t) const {return m_sigmaBar * a_St;};

		double GetStartPoint() const {return m_S0;};
	};
}
