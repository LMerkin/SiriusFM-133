#include "IRProviderConst.h"
#include <cstdio>
#include <cstdlib>
#include <iostream>

#define BUF_SIZE 1024
#define CCY_SIZE 3

namespace SiriusFM
{
	IRProvider<IRModeE::Const>::IRProvider(const char* a_file)
	{
		FILE* src = fopen(a_file, "r");
		char buf[BUF_SIZE];
		char ccy[CCY_SIZE+1];

		for(int k = 0; k < int(CcyE::N); ++k)
			m_IRs[k] = 0;

		if(a_file == nullptr) //check if a_file empty
			return;

		if(!src)
			throw std::invalid_argument("Constructor");

		while(fgets(ccy, CCY_SIZE+1, src))
		{
			fgets(buf, BUF_SIZE, src);
			m_IRs[int(Str2CcyE(ccy))] = strtod(buf+2, nullptr);//consider buf+1
		}

	}
}
