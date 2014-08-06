#include "bignum.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main(int argc, char *argv[])
{
	(void)argc;
	(void)argv;
	srand(unsigned(time(NULL)));

	BigNum bnX = 1.5f;
	float fX = bnX;
	BigNum bnY = bnX;
	BigNum bnZ = bnX+bnY;
	float fZ = bnZ;
	BigNum bnN = -bnX;
	BigNum bnN2 = -bnN;
	float fN = bnN;
	return 0;
}