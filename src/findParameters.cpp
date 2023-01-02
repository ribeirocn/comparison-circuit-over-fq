/*
* Code to adjust parameters for the Private Set Meembership primitive
*/
#include <iostream>
#include <time.h>
#include <random>

#include <helib/helib.h>
#include <helib/debugging.h>
#include <helib/Context.h>
#include <helib/polyEval.h>

#include "../../HElib/src/PrimeGenerator.h"
#include "tools.h"
#include "comparator.h"

using namespace std;
using namespace NTL;
using namespace helib;
using namespace he_cmp;


long logProduct(std::set<long> list) {
  long result = 0;
  for(long p : list ) {
    result += log(p);
  }
  return result;
}

long ctxtPrimeSize(long nBits)
{
  double bit_loss =
      -std::log1p(-1.0 / double(1L << helib::PrimeGenerator::B)) / std::log(2.0);
  // std::cerr << "*** bit_loss=" << bit_loss;

  // How many primes of size HELIB_SP_NBITS it takes to get to nBits
  double maxPsize = HELIB_SP_NBITS - bit_loss;
  // primes of length len are guaranteed to be at least (1-1/2^B)*2^len,

  long nPrimes = long(ceil(nBits / maxPsize));
  // this is sufficiently many primes

  // now we want to trim the size to avoid unnecessary overshooting
  // so we decrease targetSize, while guaranteeing that
  // nPrimes primes of length targetSize multiply out to
  // at least nBits bits.

  long targetSize = HELIB_SP_NBITS;
  while (10 * (targetSize - 1) >= 9 * HELIB_SP_NBITS &&
         (targetSize - 1) >= 30 &&
         ((targetSize - 1) - bit_loss) * nPrimes >= nBits)
    targetSize--;

  if (((targetSize - 1) - bit_loss) * nPrimes >= nBits)
    Warning(__func__ + std::string(": non-optimal targetSize"));

  return targetSize;
}

void addSpecialPrimes(std::set<long> &primes, unsigned long m, unsigned long p, unsigned long qs) {

  long phim = phi_N(m);
  long nDgts = 3;
  long nCtxtPrimes = primes.size();
  double stdev = to_double(3.2);
  if (nDgts > nCtxtPrimes)
    nDgts = nCtxtPrimes; // sanity checks
  if (nDgts <= 0)
    nDgts = 1;
  std::vector<std::set<long>> digits;
  digits.resize(nDgts); // allocate space

  if (nDgts > 1) { // we break ciphertext into a few digits when key-switching
    // NOTE: The code below assumes that all the ctxtPrimes have roughly the
    // same size

    std::set<long> remaining = primes;
    for (long dgt = 0; dgt < nDgts - 1; dgt++) {
      long digitCard = divc(remaining.size(), nDgts - dgt);
      // ceiling(#-of-remaining-primes, #-or-remaining-digits)

      for (long i : remaining) {
        digits[dgt].insert(i);
        if (digits[dgt].size() >= digitCard)
          break;
      }
      for(long i : digits[dgt]) {
        remaining.erase(i);
      } 
    }
    // The last digit has everything else
    if (empty(remaining)) { // sanity check, use one less digit
      nDgts--;
      digits.resize(nDgts);
    } else
      digits[nDgts - 1] = remaining;
  } else { // only one digit
    digits[0] = primes;
  }

  double maxDigitLog = 0.0;
  for (auto& digit : digits) {
    double size = logProduct(digit);
    if (size > maxDigitLog)
      maxDigitLog = size;
  }

  // Add special primes to the chain for the P factor of key-switching
  double nBits;

  double h = phim / 2.0;
  double log_phim = std::log(phim);
  if (log_phim < 1)
    log_phim = 1;

  if (floor(log2(m)) == ceil(log2(m))) { // is power of 2
      nBits = (maxDigitLog + std::log(stdev) +
               0.5 * std::log(12.0) + std::log(nDgts) -
               0.5 * std::log(log_phim) - 2 * std::log(p) - std::log(h)) /
              std::log(2.0);
  } else {
      nBits =
          (maxDigitLog + std::log(m) + std::log(stdev) +
           0.5 * std::log(12.0) + std::log(nDgts) - 0.5 * log_phim -
           0.5 * std::log(log_phim) - 2 * std::log(p) - std::log(h)) /
          std::log(2.0);
    
  }

    // Both of the above over-estimate nBits by a factor of
    // log2(scale). That should provide a sufficient safety margin.
    // See design document

  if (nBits < 1)
    nBits = 1;

  double bit_loss =
      -std::log1p(-1.0 / double(1L << PrimeGenerator::B)) / std::log(2.0);

  // How many primes of size HELIB_SP_NBITS it takes to get to nBits
  double maxPsize = HELIB_SP_NBITS - bit_loss;
  // primes of length len are guaranteed to be at least (1-1/2^B)*2^len,

  long nPrimes = long(ceil(nBits / maxPsize));
  // this is sufficiently many prime

  // now we want to trim the size to avoid unnecessary overshooting
  // so we decrease targetSize, while guaranteeing that
  // nPrimes primes of length targetSize multiply out to
  // at least nBits bits.

  long targetSize = HELIB_SP_NBITS;
  while ((targetSize - 1) >= 0.55 * HELIB_SP_NBITS && (targetSize - 1) >= 30 &&
         ((targetSize - 1) - bit_loss) * nPrimes >= nBits)
    targetSize--;

  if (((targetSize - 1) - bit_loss) * nPrimes >= nBits)
    Warning(__func__ + std::string(": non-optimal targetSize"));

  PrimeGenerator gen(targetSize, m);

  while (nPrimes > 0) {
    long q = gen.next();
    if (primes.count(p))
      continue;
    // nbits could equal NTL_SP_BITS or the size of one
    // of the small primes, so we have to check for duplicates here...
    // this is not the most efficient way to do this,
    // but it doesn't make sense to optimize this any further

    primes.insert(q);
    nPrimes--;
  }
}

uint calculateSecurityLevel(unsigned long p, unsigned long m, unsigned long qs) {

  long target = ctxtPrimeSize(qs);

  std::set<long> primes;
  helib::PrimeGenerator gen(target, m);
  double bitlen = 0; // how many bits we already have
  while (bitlen < qs - 0.5) {
    long q = gen.next(); // generate the next prime
    primes.insert(q);     // add it to the list
    bitlen += std::log2(q);
  }
  addSpecialPrimes(primes, m, p, qs);
  double s = to_double(3.2);
  if (floor(log2(m)) != ceil(log2(m))) { // not power of two
    s *= sqrt(m);
  }
  double log2AlphaInv = (logProduct(primes) - log(s)) / log(2.0);
  return lweEstimateSecurity(phi_N(m), log2AlphaInv, 0);
}

void adjustingParameters(unsigned long& p, unsigned long& m, unsigned long nb_primes, unsigned long d, long ss_size=-1) {
    cout << "Ajusting parameters for PSM ... " << endl;
    bool adjust_ss = ss_size < 0;
    unsigned long p_ = p;
    unsigned long m_ = m;
    // find prime p > inputed p
    for(uint i=0; i<100; i++, p++) {
      if (NTL::ProbPrime(p, 60)) {
        std::vector<long> gens;
        std::vector<long> ords;

        for(uint i=0; i<100000; i++, m++) {
          long ordP = helib::multOrd(p,m);
          if(adjust_ss) { ss_size = (p-1)>>1;};
          cout << "P: " << p <<  " M: " << m << " OrdP: " << ordP << endl;
          if(m % p != 0 && ordP >= d && ordP <  25/*&& (ordP < d + 6 || phi_N(m) >= ordP * ss_size )*/) {
            helib::findGenerators(gens,ords,m,p);
            cout << "G: " << gens.size() << endl;
            if( gens.size() < 16 /*|| gens[0] == 2*/) {
              cout << "S: " << calculateSecurityLevel(p,m,nb_primes) << endl;
              if(calculateSecurityLevel(p,m,nb_primes)>120) {
                return;
              }
            }
          }
        }
      }
    }
    cout << "Fail ajusting parameters, reusing original ones" << endl;
    p = p_;
    m = m_;
}