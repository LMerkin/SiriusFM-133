#include "IRProviderConst.h"
#include <cstdio>
#include <cstdlib>
#include <iostream>

namespace SiriusFM
{
	IRProvider<IRModeE::Const>::IRProvider(const char* a_file)
	{
    constexpr int BUF_SIZE = 64;
    constexpr int CCY_SIZE = 3;

    // Zero-out all the rates:
		for(int k = 0; k < int(CcyE::N); ++k)
			m_IRs[k] = 0;

    if (a_file == nullptr || *a_file == '\0')
      return;

		FILE* src = fopen(a_file, "r");
		if(a_file == nullptr)
      // File is non-existent:
      throw std::runtime_error("Cannot open file");

		char buf[BUF_SIZE];
		char ccy[CCY_SIZE + 1] = "XXX";

		while(fgets(buf, BUF_SIZE, src) != nullptr)
		{
      if (*buf == '\0' || *buf == '\n' || *buf == '#')
        continue;

      strncpy(ccy, buf, CCY_SIZE);
			m_IRs[int(Str2CcyE(ccy))] = atof(buf + CCY_SIZE + 1);
		}
    fclose(src);
	}
}
