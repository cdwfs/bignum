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
		m_intLength  = 1; // TEMP;
		m_fracLength = 3; // TEMP;
	}
	BigNum(const BigNum &rhs)
	{
		*this = rhs;
	}
	BigNum& operator=(const BigNum& rhs)
	{
		if (this != &rhs)
		{
			m_sign = rhs.m_sign;
			assert(rhs.m_intLength == 1 && rhs.m_fracLength == 3); // TEMP
			m_intLength  = rhs.m_intLength;
			m_fracLength = rhs.m_fracLength;
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
		m_sign = fbits.s;
		m_intLength  = 1; // TEMP
		m_fracLength = 3; // TEMP
		int32_t exponent = (int32_t)fbits.e - 127;
		exponent += 1; // Correct for implicit leading 1
		if (exponent < -3*32)
		{
			// TODO: too conservative; ignores denormalized values, but we can't represent them anyway in 32.96
			return;
		}
		if (exponent > 32)
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
		m_value[iDW1] = converter.hi32;
		if (iDW2 < 4)
		{
			m_value[iDW2] = converter.lo32;
		}
	}
	operator float() const
	{
		// Find total number of leading zeroes
		const int32_t kMaxExponent = 32; // largest possible for 32.96
		int32_t exponent = kMaxExponent;
		int iDW;
		for(iDW=0; iDW<4; ++iDW)
		{
			int32_t lz = (int32_t)countLeadingZeroes(m_value[iDW]);
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
		converter.hi32  = m_value[iDW1];
		converter.lo32 = (iDW2<4) ? m_value[iDW2] : 0;
		converter.asQword >>= ((64-23) - countLeadingZeroes(converter.hi32) - 1);
		fbits.m = converter.mantissa;
		fbits.e = exponent + 127; // 32-bit float bias
		fbits.s = m_sign;
		return fbits.f;
	}

	BigNum& operator+=(const BigNum &rhs)
	{
		if (m_sign == rhs.m_sign) // signs match: compute lhs+rhs, keep lhs's sign
		{
			uint8_t carry = 0;
			for(int iDW=3; iDW>=0; --iDW)
			{
				carry = _addcarry_u32(carry, m_value[iDW], rhs.m_value[iDW], m_value+iDW);
			}
		}
		else if (abs(rhs) < abs(*this)) // signs differ, and |lhs| > |rhs|: compute lhs-rhs and keep lhs's sign
		{
			uint8_t borrow = 0;
			for(int iDW=3; iDW>=0; --iDW)
			{
				borrow = _subborrow_u32(borrow, m_value[iDW], rhs.m_value[iDW], m_value+iDW);
			}
		}
		else // signs differ and |lhs| <= |rhs|: compute rhs-lhs and flip lhs's sign
		{
			uint8_t borrow = 0;
			for(int iDW=3; iDW>=0; --iDW)
			{
				borrow = _subborrow_u32(borrow, rhs.m_value[iDW], m_value[iDW], m_value+iDW);
			}
			m_sign = 1-m_sign;
		}
		return *this;
	}

	BigNum& operator*=(const BigNum &rhs)
	{
		uint32_t product[8] = {0}; // TODO: should be double the size of a regular bignum

		// Implements basic "schoolbook" multiplication. More efficient algorithms exist, but they're
		// not faster until m_length >10ish
		for(int iR=3; iR>=0; --iR)
		{
			uint8_t carryT = 0, carryO = 0;
			uint32_t hiPrev = 0;
			for(int iL=3; iL>=0; --iL)
			{
				uint32_t temp, lo, hi;
				lo = _mulx_u32(m_value[iL], rhs.m_value[iR], &hi); // [hi,lo] = lhs[iL] * rhs[iR]
				carryT = _addcarry_u32(carryT, lo, hiPrev, &temp); // [temp, carryT] = lo + hiPrev + carryT
				hiPrev = hi;
				carryO = _addcarry_u32(carryO, temp, product[1+iL+iR], product+(1+iL+iR)); // [prod[1+iL+iR], carryO] += temp+carryO
			}
			carryT = _addcarry_u32(carryT, hiPrev, product[iR], product+iR); // prod[iR] += (hiPrev+carryT)
			assert(carryT == 0);
			carryO = _addcarry_u32(carryO,      0, product[iR], product+iR); // prod[iR] += carryO
			assert(carryO == 0);
		}
		m_sign ^= rhs.m_sign;
		m_value[0] = product[1];
		m_value[1] = product[2];
		m_value[2] = product[3];
		m_value[3] = product[4];
		return *this;
	}

	// Comparison
	bool operator==(const BigNum &rhs) const
	{
		// TODO: handling different lengths is tricky.
		assert(m_intLength == rhs.m_intLength && m_fracLength == rhs.m_fracLength);
		bool isZero = true;
		for(uint32_t iDW=0; iDW<m_intLength+m_fracLength; ++iDW)
		{
			if (m_value[iDW] != rhs.m_value[iDW])
				return false; // value bits mismatch
			isZero = isZero && (m_value[iDW] == 0);
		}
		if (isZero)
			return true; // value bits match, and all are zero. Ignore the sign bit.
		return m_sign == rhs.m_sign; // magnitudes are equal and non-zero. Do the sign bits match?
	}
	bool operator<(const BigNum &rhs) const
	{
		// TODO: handling different lengths is tricky.
		assert(m_intLength == rhs.m_intLength && m_fracLength == rhs.m_fracLength);
		bool isZero = 0;
		for(uint32_t iDW=0; iDW<m_intLength+m_fracLength; ++iDW)
		{
			if      (m_value[iDW] < rhs.m_value[iDW]) return (rhs.m_sign == 0);
			else if (m_value[iDW] > rhs.m_value[iDW]) return (m_sign != 0);
			isZero = isZero && ((m_value[iDW] | rhs.m_value[iDW]) == 0);
		}
		// |lhs| == |rhs|
		// lhs < rhs iff |lhs| != 0 and lhs is negative and rhs is not
		return !isZero && (m_sign > rhs.m_sign);
	}

	const BigNum operator+(void) const
	{
		return *this;
	}
	const BigNum operator-(void) const
	{
		BigNum out = *this;
		out.m_sign = ~m_sign;
		return out;
	}

	uint32_t m_intLength : 15;  // in DWORDS
	uint32_t m_fracLength : 16; // in DWORDS
	uint32_t m_sign : 1;        // 0: positive, 1: negative
	uint32_t m_value[4]; // precision is implicitly 32.96 for now

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
inline bool operator!=(const BigNum &lhs, const BigNum &rhs) { return !(lhs == rhs); }
inline bool operator> (const BigNum &lhs, const BigNum &rhs) { return  (rhs  < lhs); }
inline bool operator<=(const BigNum &lhs, const BigNum &rhs) { return !(lhs  > rhs); }
inline bool operator>=(const BigNum &lhs, const BigNum &rhs) { return !(lhs  < rhs); }
inline BigNum operator+(BigNum lhs, const BigNum &rhs) { lhs +=  rhs; return lhs; }
inline BigNum operator-(BigNum lhs, const BigNum &rhs) { lhs += -rhs; return lhs; }
inline BigNum operator*(BigNum lhs, const BigNum &rhs) { lhs *=  rhs; return lhs; }
inline BigNum abs(BigNum bn) { bn.m_sign = 0; return bn; }
#pragma warning(default:4201) // nameless struct/union
