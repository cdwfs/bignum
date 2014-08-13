#include "bignum.h"

#include <algorithm>
#include <assert.h>
#include <intrin.h>
#include <limits>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#pragma warning(disable:4201) // nameless struct/union

static inline uint32_t countLeadingZeroes(const BigNum bn)
{
	unsigned long pos    = 0;
	unsigned long offset = 0;
	for(int iDW=0; iDW<(sizeof(bn)/sizeof(uint32_t)); ++iDW)
	{
		uint8_t result = _BitScanReverse(&pos, bn.m_value[iDW]);
		if (result)
			return offset + (31-pos);
		offset += 32;
	}
	return offset;
}


int main(int argc, char *argv[])
{
	(void)argc;
	(void)argv;
	srand(unsigned(time(NULL)));

	//
	// Test float -> BigNum -> float round-trip conversion
	//
	union {
		float f;
		uint32_t i;
		struct {
			uint32_t m : 23;
			uint32_t e:   8;
			uint32_t s:   1;
		};
	} fbitsIn, fbitsOut;

	BigNum bnMin, bnMax;
	memset(&bnMin, 0x00, sizeof(bnMin));
	memset(&bnMax, 0xFF, sizeof(bnMax));
	bnMin.m_value[3] = 1;
	bnMax.m_value[0] = 0x7FFFFFFF;
	float fMin = bnMin; // 0x0F800000
	float fMax = bnMax; // 0x4EFFFFFF

	unsigned int fpcOld = 0, fpcNew = 0;
	_controlfp_s(&fpcOld, 0,0); // retrieve current FP control bits
	_controlfp_s(&fpcNew, _RC_CHOP, _MCW_RC); // Set rounding mode to chop (round towards zero); much easier to emulate in software than round-to-nearest.

	uint32_t numTests = 0;
	uint32_t numErrors = 0;
	for(uint32_t iF=0; iF<0xFFFFFFFF; ++iF)
	{
		if ((iF & 0xFFFFF) == 0)
		{
			printf("0x%08X (%u tests, %u errors)\n", iF, numTests, numErrors);
		}
		fbitsIn.i = iF;
		if (fbitsIn.f != fbitsIn.f  || // Can't convert NaN
			fabsf(fbitsIn.f) < fMin || fabsf(fbitsIn.f) > fMax)
		{
			continue;
		}
		BigNum bn = fbitsIn.f;
		fbitsOut.f = bn;
		// When converting an extremely small float to BigNum, we're throwing away
		// all but the top bit(s) of the mantissa. Converting back is thus necessarily quite lossy.
		// if numLostBits < 24, we need to mask bits off
		BigNum absBN = (fbitsIn.f < 0) ? -bn : bn;
		uint32_t highestSetBitIndex = 8*sizeof(absBN) - countLeadingZeroes(absBN);
		int32_t numLostBits = std::max(0, 24-int32_t(highestSetBitIndex));
		uint32_t lostBitsMask = (1 << numLostBits)-1;
		fbitsIn.i &= ~lostBitsMask;
		// an input of -0.0 will convert to +0.0; make sure that counts as a positive match
		if (fbitsIn.f == -0.0f)
			fbitsIn.f = +0.0f;
		if ( abs(int(fbitsOut.i-fbitsIn.i)) > 0) // differ by more than one bit
		{
			fprintf(stderr, "Conversion error: %.9e (0x%08X) -> BigNum -> %.9e (0x%08X)\n",
				fbitsIn.f, fbitsIn.i, fbitsOut.f, fbitsOut.i);
			++numErrors;
			__debugbreak();
		}
		++numTests;
	}

	_controlfp_s(&fpcNew, fpcOld, _MCW_RC); // Restore original rounding mode
}
