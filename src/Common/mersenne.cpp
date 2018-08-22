/************************** MERSENNE.CPP ******************** AgF 2001-10-18 *
*  Random Number generator 'Mersenne Twister'                                *
*                                                                            *
*  This random number generator is described in the article by               *
*  M. Matsumoto & T. Nishimura, in:                                          *
*  ACM Transactions on Modeling and Computer Simulation,                     *
*  vol. 8, no. 1, 1998, pp. 3-30.                                            *
*  Details on the initialization scheme can be found at                      *
*  http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html                  *
*                                                                            *
*  Experts consider this an excellent random number generator.               *
*                                                                            *
*  © 2001 - 2004 A. Fog.                                                     *
*  GNU General Public License www.gnu.org/copyleft/gpl.html                  *
*****************************************************************************/

/*------------------------------------------------------------------------------

Modified for RangeShifter v2.0

Last updated: 5 December 2016 by Steve Palmer

------------------------------------------------------------------------------*/

#include "randomc.h"

//---------------------------------------------------------------------------
// RS version control added by SCFP 23/8/16
#include "Version.h"
#if RSDEBUG
#include "Parameters.h"
extern paramSim *paramsSim;
//#include <string>
//#include <stdio.h>
#include <fstream>
//#include <sstream>
//#include <iostream>
//#include <io.h>
//#include <iomanip>
#include <stdlib.h>
using namespace std;
//extern ofstream DEBUGLOG;
int limitdebuglines = 0;
#endif
//---------------------------------------------------------------------------

#if RSDEBUG
// RS random initialisation log added by SCFP 29/8/16
ofstream RANDOMLOG;
#endif

void TRandomMersenne::RandomInit(uint32 seed) { // re-seed generator

#if RSDEBUG

// RS random initialisation log added by SCFP 25/8/16
string name = paramsSim->getDir(2) + "RandomLog.txt";
RANDOMLOG.open(name.c_str());
RANDOMLOG << "TRandomMersenne::RandomInit(): seed=" << seed << endl;

// NOTE: hex value 3FF00000 is the same as decimal integer 1072693248

#endif
	mt[0]= seed;
	for (mti=1; mti < MERS_N; mti++) {
		mt[mti] = (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti);}

	// detect computer architecture
	union {double f; uint32 i[2];} convert;
#if RSDEBUG
RANDOMLOG << "TRandomMersenne::RandomInit(): convert.f=" << convert.f
	<< " convert.i[0]=" << convert.i[0] << " convert.i[1]=" << convert.i[1]
	<< endl;
#endif
	convert.f = 1.0;
#if RSDEBUG
RANDOMLOG << "TRandomMersenne::RandomInit(): convert.f=" << convert.f
	<< " convert.i[0]=" << convert.i[0] << " convert.i[1]=" << convert.i[1]
	<< endl;
#endif
	// Note: Old versions of the Gnu g++ compiler may make an error here,
	// compile with the option  -fenum-int-equiv  to fix the problem
	if (convert.i[1] == 0x3FF00000) {
		Architecture = LITTLE_ENDIAN1;
#if RSDEBUG
RANDOMLOG << "TRandomMersenne::RandomInit(): Architecture = LITTLE_ENDIAN1"
	<< endl;
#endif
	}
	else {
		if (convert.i[0] == 0x3FF00000) {
			Architecture = BIG_ENDIAN1;
#if RSDEBUG
RANDOMLOG << "TRandomMersenne::RandomInit(): Architecture = BIG_ENDIAN1"
	<< endl;
#endif
		}
		else {
			Architecture = NONIEEE;
#if RSDEBUG
RANDOMLOG << "TRandomMersenne::RandomInit(): Architecture = NONIEEE"
	<< endl;
#endif
		}
	}
#if RSDEBUG
RANDOMLOG.close(); RANDOMLOG.clear();
#endif
}

  
void TRandomMersenne::RandomInitByArray(uint32 seeds[], int length) {
  // seed by more than 32 bits
  int i, j, k;
  RandomInit(19650218UL);
  if (length <= 0) return;
  i = 1;  j = 0;
  k = (MERS_N > length ? MERS_N : length);
  for (; k; k--) {
    mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL)) + seeds[j] + j;
    i++; j++;
    if (i >= MERS_N) {mt[0] = mt[MERS_N-1]; i=1;}
    if (j >= length) j=0;}
  for (k = MERS_N-1; k; k--) {
    mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL)) - i;
    if (++i >= MERS_N) {mt[0] = mt[MERS_N-1]; i=1;}}
  mt[0] = 0x80000000UL;} // MSB is 1; assuring non-zero initial array

  
uint32 TRandomMersenne::BRandom() {
  // generate 32 random bits
  uint32 y;

  if (mti >= MERS_N) {
    // generate MERS_N words at one time
    const uint32 LOWER_MASK = (1LU << MERS_R) - 1; // lower MERS_R bits
    const uint32 UPPER_MASK = -1L  << MERS_R;      // upper (32 - MERS_R) bits
    static const uint32 mag01[2] = {0, MERS_A};
    
    int kk;
    for (kk=0; kk < MERS_N-MERS_M; kk++) {    
      y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
      mt[kk] = mt[kk+MERS_M] ^ (y >> 1) ^ mag01[y & 1];}

    for (; kk < MERS_N-1; kk++) {    
      y = (mt[kk] & UPPER_MASK) | (mt[kk+1] & LOWER_MASK);
      mt[kk] = mt[kk+(MERS_M-MERS_N)] ^ (y >> 1) ^ mag01[y & 1];}      

    y = (mt[MERS_N-1] & UPPER_MASK) | (mt[0] & LOWER_MASK);
    mt[MERS_N-1] = mt[MERS_M-1] ^ (y >> 1) ^ mag01[y & 1];
    mti = 0;}

  y = mt[mti++];

  // Tempering (May be omitted):
  y ^=  y >> MERS_U;
  y ^= (y << MERS_S) & MERS_B;
  y ^= (y << MERS_T) & MERS_C;
  y ^=  y >> MERS_L;
  return y;}

  
double TRandomMersenne::Random() {
  // output random double number in the interval 0 <= x < 1
	union {double f; uint32 i[2];} convert;
	uint32 r = BRandom(); // get 32 random bits
  // The fastest way to convert random bits to floating point is as follows:
  // Set the binary exponent of a floating point number to 1+bias and set
  // the mantissa to random bits. This will give a random number in the 
  // interval [1,2). Then subtract 1.0 to get a random number in the interval
  // [0,1). This procedure requires that we know how floating point numbers
  // are stored. The storing method is tested in function RandomInit and saved 
  // in the variable Architecture. The following switch statement can be
	// omitted if the architecture is known. (A PC running Windows or Linux uses
  // LITTLE_ENDIAN1 architecture):
	switch (Architecture) {
	case LITTLE_ENDIAN1:
#if RSDEBUG
//if (limitdebuglines < 10) {
//	DEBUGLOG << "TRandomMersenne::Random(): case LITTLE_ENDIAN1 "
//		<< " r=" << r << endl;
//}
#endif
		convert.i[0] =  r << 20;
#if RSDEBUG
//if (limitdebuglines < 10) {
//	DEBUGLOG << "TRandomMersenne::Random(): " << " r=" << r
//		<< " convert.i[0]=" << convert.i[0] << " convert.f=" << convert.f << endl;
//}
#endif
		convert.i[1] = (r >> 12) | 0x3FF00000;
#if RSDEBUG
//if (limitdebuglines < 10) {
//	DEBUGLOG << "TRandomMersenne::Random(): " << " r=" << r
//		<< " convert.i[1]=" << convert.i[1] << " convert.f=" << convert.f << endl;
//	limitdebuglines++;
//}
#endif
#if RSDEBUG
//RANDOMLOG << convert.f - 1.0 << endl;
#endif
		return convert.f - 1.0;
	case BIG_ENDIAN1:
#if RSDEBUG
//if (limitdebuglines < 10) {
//	DEBUGLOG << "TRandomMersenne::Random(): case BIG_ENDIAN1 " << endl;
//	limitdebuglines++;
//}
#endif
		convert.i[1] =  r << 20;
		convert.i[0] = (r >> 12) | 0x3FF00000;
		return convert.f - 1.0;
	case NONIEEE:
#if RSDEBUG
//if (limitdebuglines < 1000) {
//	DEBUGLOG << "TRandomMersenne::Random(): case NONIEEE "
//		<< " r=" << r << endl;
//}
#endif
		union {double fff; uint32 iii;} linuxconvert;
		linuxconvert.iii = (r >> 12) | 0x3FF0000000000000;
#if RSDEBUG
//if (limitdebuglines < 1000) {
//	DEBUGLOG << "TRandomMersenne::Random(): " << " r=" << r
//		<< " linuxconvert.iii=" << linuxconvert.iii
//		<< " linuxconvert.fff=" << linuxconvert.fff << endl;
//	limitdebuglines++;
//}
#endif
#if RSDEBUG
//RANDOMLOG << linuxconvert.fff - 1.0 << endl;
#endif
		return linuxconvert.fff - 1.0;
	default:
  ;} 
  // This somewhat slower method works for all architectures, including 
  // non-IEEE floating point representation:
  return (double)r * (1./((double)(uint32)(-1L)+1.));}

  
int TRandomMersenne::IRandom(int min, int max) {
  // output random integer in the interval min <= x <= max
  int r; 
  r = int((max - min + 1) * Random()) + min; // multiply interval with random and truncate
  if (r > max) r = max;
  if (max < min) return 0x80000000;
  return r;}

