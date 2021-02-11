#pragma once
#include<cmath>

namespace SiriusFM
{
                
	class DiffusionOU
	{
		double const m_theta;
		double const m_kappa;
		double const m_sigma;

	public:
		DiffusionOU(double m, double s, double t): 	m_theta(m),
													m_kappa(s), 
													m_sigma(t)
		{
			if(m_sigma <= 0)
			{
			}
		};

		double mu(double S_t, double t) {return m_kappa*(m_theta -  S_t);};
		double sigma(double S_t = 0, double t) {return m_sigma;};
	};
}
