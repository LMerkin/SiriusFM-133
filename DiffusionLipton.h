#pragma once
namespace SiriusFM
{
	class DiffusionLipton
	{
		double const m_mu;
		double const m_sigma_0;
		double const m_sigma_1;
		double const m_sigma_2;

	public:

		DiffusionLipton(double m, double s0, double s1, double s2):
			m_mu(m),
			m_sigma_0(s0),
			m_sigma_1(s1),
			m_sigma_2(s2)
		{
			if(m_sigma_1 * m_sigma_1 >= 4 * m_sigma_0 * m_sigma_2)
			{
			}
		};

		double mu(double s, double t){return m_mu * s;}
		double sigma(double s, double t)
		{
			return m_sigma_0 + m_sigma_1 * s + m_sigma_2 * s * s;
		}
	};
}
