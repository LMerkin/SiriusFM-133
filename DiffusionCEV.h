#pragma once
#include<cmath>

namespace SiriusFM
{
                
	class DiffusionCEV
	{
		double const m_muBar;
		double const m_sigmaBar;
		double const m_beta;

	public:
		DiffusionCEV(double m, double s, double b): m_muBar(m), 
											  		m_sigmaBar(s),
													m_beta(b)
		{
			if(m_sigmaBar <= 0 || m_beta <= 0)
			{
			}
		};

		double mu(double S_t, double t) {return m_muBar * S_t;};
		double sigma(double S_t, double t) 
		{
			return m_sigmaBar * pow(S_t, m_beta);
		};
	};
}
