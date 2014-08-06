#include <assert.h>
#include <intrin.h>
#include <limits>
#include <stdint.h>

#pragma warning(disable:4201) // nameless struct/union
class BigNum
{
public:
	BigNum(void)
	{
	}
	BigNum(const BigNum &rhs)
	{
		m_value[0] = rhs.m_value[0];
		m_value[1] = rhs.m_value[1];
		m_value[2] = rhs.m_value[2];
		m_value[3] = rhs.m_value[3];
	}
	BigNum& operator=(const BigNum& rhs)
	{
		if (this != &rhs)
		{
			m_value[0] = rhs.m_value[0];
			m_value[1] = rhs.m_value[1];
			m_value[2] = rhs.m_value[2];
			m_value[3] = rhs.m_value[3];
		}
		return *this;
	}

	//
	// Conversion operators
	//
	BigNum(float f)
	{
		m_value[0] = m_value[1] = m_value[2] = m_value[3] = 0;
		SingleBits fbits;
		fbits.f = f;
		int32_t exponent = (int32_t)fbits.e - 127;
		exponent += 1; // Correct for implicit leading 1
		if (exponent < -3*32)
		{
			// TODO: too conservative; ignores denormalized values, but we can't represent them anyway in 32.96
			return;
		}
		if (exponent >= 32)
		{
			// Too large to represent!
			assert(0);
			m_value[0] = 0x7FFFFFFF;
			m_value[1] = 0xFFFFFFFF;
			m_value[2] = 0xFFFFFFFF;
			m_value[3] = 0xFFFFFFFF;
			return;
		}
		Mantissa64 converter;
		converter.unused   = 1; // implicit 1.mantissa
		converter.mantissa = fbits.m;
		// exponent is [31..-96]
		const int32_t kMaxExponent = 32;
		const uint32_t lzTotal = kMaxExponent-exponent;
		int32_t iDW1 = lzTotal / 32;
		int32_t iDW2 = iDW1+1;
		converter.asQword <<= ((64-24) - (lzTotal % 32));
		m_value[iDW1] = int32_t(converter.hi32);
		if (iDW2 < 4)
		{
			m_value[iDW2] = int32_t(converter.lo32);
		}
		if (f < 0)
		{
			*this = -*this;
		}
	}
	operator float() const
	{
		BigNum av = *this;
		if (m_value[0] >> 31)
		{
			// value is negative; convert local copy to positive
			av = -av;
		}
		// Find total number of leading zeroes
		const int32_t kMaxExponent = 32; // largest possible for 32.96
		int32_t exponent = kMaxExponent;
		int iDW;
		for(iDW=0; iDW<4; ++iDW)
		{
			int32_t lz = (int32_t)countLeadingZeroes(av.m_value[iDW]);
			exponent -= lz;
			if (lz < 32)
			{
				break; // found a one!
			}
		}
		if (exponent == -3*32)
		{
			return 0; // all zeroes
		}
		else if (exponent <= -126)
		{
			assert(0); // TODO: denormalized; possibly zero.
			return 0;
		}
		else if (exponent > 127)
		{
			assert(0); // TODO: too large for destination float; clamp to infinity
			return std::numeric_limits<float>::infinity();
		}
		exponent -= 1; // skip implicit leading 1

		SingleBits fbits;
		int32_t iDW1 = iDW;
		int32_t iDW2 = iDW1+1;
		Mantissa64 converter;
		converter.hi32  = av.m_value[iDW1];
		converter.lo32 = (iDW2<4) ? av.m_value[iDW2] : 0;
		converter.asQword >>= ((64-23) - countLeadingZeroes(converter.hi32) - 1);
		fbits.m = converter.mantissa;
		fbits.e = exponent + 127; // 32-bit float bias
		fbits.s = (m_value[0] >> 31) & 0x1; // use original value for sign bit!
		return fbits.f;
	}

	BigNum& operator+=(const BigNum &rhs)
	{
		uint64_t carry = 0;
		for(int iDW=3; iDW>=0; --iDW)
		{
			uint64_t result = carry + uint64_t(uint32_t(this->m_value[iDW])) + uint64_t(uint32_t(rhs.m_value[iDW]));
			this->m_value[iDW] = int32_t(result); // implicit & 0xFFFFFFFF
			carry = result >> 32;
			assert(carry <= 1); // carry bit must be either 0 or 1
		}
		assert(carry == 0); // If there's carry in the first DW, we have nowhere to put it
		return *this;
	}

	const BigNum operator+(void) const
	{
		return *this;
	}
	const BigNum operator-(void) const
	{
		BigNum out;
		uint64_t carry = 1;
		for(int iDW=3; iDW>=0; --iDW)
		{
			uint64_t result = carry + uint64_t(uint32_t(~m_value[iDW]));
			out.m_value[iDW] = int32_t(result); // implicit & 0xFFFFFFFF
			carry = result >> 32;
			assert(carry <= 1);
		}
		return out;
	}

	int32_t m_value[4]; // precision is implicitly 32.96 for now

private:
	typedef union { float f;  struct { uint32_t m:23; uint32_t e: 8; uint32_t s:1; }; } SingleBits;
	typedef union { double d; struct { uint64_t m:52; uint64_t e:11; uint64_t s:1; }; } DoubleBits;
	typedef union
	{
		struct
		{
			uint64_t lo32 : 32;
			uint64_t hi32  : 32;
		};
		uint64_t asQword;
		struct
		{
			uint64_t mantissa : 23;
			uint64_t unused   : 41;
		};
	} Mantissa64;

	static inline uint32_t countLeadingZeroes(const uint32_t in)
	{
		unsigned long pos = 0;
		uint8_t result = _BitScanReverse(&pos, in);
		return result ? (31-pos) : 32;
	}
};
inline BigNum operator+(BigNum lhs, const BigNum &rhs) { lhs += rhs; return lhs; }
#pragma warning(default:4201) // nameless struct/union
