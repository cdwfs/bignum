#include "bignum.h"

#include <algorithm>
#include <assert.h>
#include <intrin.h>
#include <limits>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <random>

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


typedef union
{
	float f;
	uint32_t i;
	struct {
		uint32_t m : 23;
		uint32_t e:   8;
		uint32_t s:   1;
	};
} SingleBits;

int main(int argc, char *argv[])
{
	(void)argc;
	(void)argv;
	std::mt19937 rng;

	unsigned int fpcOld = 0, fpcNew = 0;
	_controlfp_s(&fpcOld, 0,0); // retrieve current FP control bits
	_controlfp_s(&fpcNew, _RC_CHOP, _MCW_RC); // Set rounding mode to chop (round towards zero); much easier to emulate in software than round-to-nearest.

	// Compute the largest and smallest float that can be represented in a BigNum
	BigNum bnMin = 0.0f, bnMinExact = 0.0f, bnMax = 0.0f;
	memset(&bnMin.m_value,      0x00, (bnMin.m_intLength+bnMin.m_fracLength)*sizeof(uint32_t));
	memset(&bnMinExact.m_value, 0x00, (bnMin.m_intLength+bnMin.m_fracLength)*sizeof(uint32_t));
	memset(&bnMax.m_value,      0xFF, (bnMin.m_intLength+bnMin.m_fracLength)*sizeof(uint32_t));
	bnMin.m_value[3] = 1;
	bnMinExact.m_value[3] = 1 << 23;
	const float fMin = bnMin; // 0x0F800000
	const float fMinExact = bnMinExact;
	const float fMax = bnMax; // 0x4F7FFFFF

	//
	// Test BigNum arithmetic
	//
	uint32_t numArithTests  = 0;
	uint32_t numArithErrors = 0;
	for(uint32_t iF=0; iF<0xFFFFFFFF; ++iF)
	{
		if ((iF & 0xFFFFF) == 0)
		{
			printf("0x%08X (%u tests, %u errors)\n", iF, numArithTests, numArithErrors);
		}
		SingleBits fx, fy;
		do
		{
			fx.i = rng();
		}
		while (fx.f != fx.f || (fx.f) < fMinExact || fabsf(fx.f) > fMax);
		do
		{
			fy.i = rng();
		}
		while (fy.f != fy.f || fabsf(fy.f) < fMinExact || fabsf(fy.f) > fMax);
		//fx.i = 0x1c1fe652;
		//fy.i = 0xcc50aa59;
		BigNum bnx = fx.f, bny = fy.f;
		++numArithTests;

		const float  sumF = fx.f + fy.f;
		const float diffF = fx.f - fy.f;
		const float prodF = fx.f * fy.f;
		const bool    ltF = fx.f < fy.f;
		const bool    eqF = fx.f == fy.f;

		// Test addition
		if (fabsf(sumF) <= fMax)
		{
			const BigNum bnSum = bnx + bny;
			const float sumB = bnSum;
			if (fabsf(sumF) < fMinExact)
			{
				if (fabsf(sumF) < fMin)
				{
					assert(sumB == 0); // small results must be clamped properly to zero
				}
				// Anything less than fMinExact isn't guaranteed to match exactly anyway,
				// so I'll ignore them for now.
			}
			else
			{
				SingleBits bitsB, bitsF;
				bitsB.f = sumB;
				bitsF.f = sumF;
				if ( abs(int(bitsB.i-bitsF.i)) > 0) // differ by more than one bit
				{
					fprintf(stderr, "Math error: %.9e (0x%08X) + %.9e (0x%08X) = %.9e (0x%08X) vs. %.9e (0x%08X)\n",
						fx.f, fx.i, fy.f, fy.i, bitsB.f, bitsB.i, bitsF.f, bitsF.i);
					++numArithErrors;
					__debugbreak();
				}
			}
		}

		// Test subtraction
		if (fabsf(diffF) <= fMax)
		{
			const BigNum bnDiff = bnx - bny;
			const float diffB = bnDiff;
			if (fabsf(diffF) < fMinExact)
			{
				if (fabsf(diffF) < fMin)
				{
					assert(diffB == 0); // small results must be clamped properly to zero
				}
				// Anything less than fMinExact isn't guaranteed to match exactly anyway,
				// so I'll ignore them for now.
			}
			else
			{
				SingleBits bitsB, bitsF;
				bitsB.f = diffB;
				bitsF.f = diffF;
				if ( abs(int(bitsB.i-bitsF.i)) > 0) // differ by more than one bit
				{
					fprintf(stderr, "Math error: %.9e (0x%08X) - %.9e (0x%08X) = %.9e (0x%08X) vs. %.9e (0x%08X)\n",
						fx.f, fx.i, fy.f, fy.i, bitsB.f, bitsB.i, bitsF.f, bitsF.i);
					++numArithErrors;
					__debugbreak();
				}
			}
		}

		// Test multiplication
		if (fabsf(prodF) <= fMax)
		{
			const BigNum bnProd = bnx * bny;
			const float prodB = bnProd;
			if (fabsf(prodF) < fMinExact)
			{
				if (fabsf(prodF) < fMin)
				{
					assert(prodB == 0); // small results must be clamped properly to zero
				}
				// Anything less than fMinExact isn't guaranteed to match exactly anyway,
				// so I'll ignore them for now.
			}
			else
			{
				SingleBits bitsB, bitsF;
				bitsB.f = prodB;
				bitsF.f = prodF;
				if ( abs(int(bitsB.i-bitsF.i)) > 0) // differ by more than one bit
				{
					fprintf(stderr, "Math error: %.9e (0x%08X) * %.9e (0x%08X) = %.9e (0x%08X) vs. %.9e (0x%08X)\n",
						fx.f, fx.i, fy.f, fy.i, bitsB.f, bitsB.i, bitsF.f, bitsF.i);
					++numArithErrors;
					__debugbreak();
				}
			}
		}

		// Test comparisons
		const bool    ltB = bnx  < bny;
		const bool    eqB = bnx == bny;
		if (ltB != ltF)
		{
			fprintf(stderr, "Comparison error: %.9e (0x%08X) < %.9e (0x%08X)\n",
				fx.f, fx.i, fy.f, fy.i);
			++numArithErrors;
			__debugbreak();
		}
		if (eqB != eqF)
		{
			fprintf(stderr, "Comparison error: %.9e (0x%08X) == %.9e (0x%08X)\n",
				fx.f, fx.i, fy.f, fy.i);
			++numArithErrors;
			__debugbreak();
		}
	}

	//
	// Test float -> BigNum -> float round-trip conversion
	//
	SingleBits fbitsIn, fbitsOut;

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
		uint32_t highestSetBitIndex = (bn.m_intLength+bn.m_fracLength)*sizeof(uint32_t)*8 - countLeadingZeroes(bn);
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
