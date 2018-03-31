#include <time.h>

#include "rng.h"
#include "SFMT-src-1.5.1/SFMT.h"

long SEED_LCG = 1, SEED_PM = 1;
sfmt_t sfmt;

void seed_lcg(long seed) {
  SEED_LCG = seed;
}

void seed_pm(long seed) {
  if (seed != 0) {
    SEED_PM = seed;
  } else {
    // To avoid seeding PM with 0
    SEED_PM = time(NULL);
  }
}

void seed_mt(int seed) {
  sfmt_init_gen_rand(&sfmt, seed);
}

// The standard lcg with parameters
// that just produce maximum period
float rand_lcg() {
  long m = 113829760;
  long a = 3557181;
  long c = 32035513;

  SEED_LCG = (SEED_LCG*a + c) % m;
  return (float) SEED_LCG / m;
}

// Park-Miller generator w/ Schrage trick
// to avoid 32-bit overflow
float rand_pm() {
  long m = 2147483647;
  long a = 16807;
  long q, r, high, low;

  q = m / a;
  r = m % a;

  high = SEED_PM / q;
  low = SEED_PM % q;

  SEED_PM = a*low - r*high;
  if (SEED_PM < 0) {
    SEED_PM += m;
  }

  return (float) SEED_PM / m;
}

// Mersenne twister
float rand_mt() {
  return (float) sfmt_genrand_uint32(&sfmt) / 4294967296;
}

// Mersenne twister integer
unsigned int rand_mt_int() {
  return sfmt_genrand_uint32(&sfmt);
}
