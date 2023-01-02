/* Copyright (C) 2019 IBM Corp.
 * This program is Licensed under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance
 * with the License. You may obtain a copy of the License at
 *   http://www.apache.org/licenses/LICENSE-2.0
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License. See accompanying LICENSE file.
 */

/*
* Extended in 2022 from https://github.com/iliailia/comparison-circuit-over-fq
* with a Private Set Membership alternative
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



// the main function that takes 7 arguments (type in Terminal: ./comparison_circuit argv[1] argv[2] argv[3] argv[4] argv[5] argv[6] argv[7] argv[8])
// argv[1] - circuit type (U, B, T or P)
// argv[2] - the plaintext modulus
// argv[3] - the dimension of a vector space over a finite field
// argv[4] - the order of the cyclotomic ring
// argv[5] - the bitsize of the ciphertext modulus in ciphertexts (HElib increases it to fit the moduli chain). The modulus used for public-key generation
// argv[6] - the length of vectors to be compared
// argv[7] - the number of experiment repetitions
// argv[8] - print debug info (y/n)

// Running examples from table 2, Section A of [Ribeiro23]
// PSM tests
// P 131 1 25743 260 1 10 y
// P 1031 1 24247 400 1 10 y
// P 2053 1 35443 440 1 10 y
// P 8209 1 39283 550 1 10 y
// P 65537 1 65536 730 1 10  32768
// Univariat Tests for comparasion
// U 131 1 25743 260 1 10	y
// U 1031 1 24247 400 1 10 y
// U 2053 1 35443 450 1 10 y
// U 8209 1 39283 560 1 10 y
// U 65537 1 65536 730 1 1 y


// extern declaration
void adjustingParameters(unsigned long& p, unsigned long& m, unsigned long nb_primes, unsigned long d, long ss_size=-11);

int main(int argc, char *argv[]) {
  if(argc < 9) {
    throw invalid_argument("There should be exactly 8 arguments\n");
  }


  CircuitType type = UNI;
  if (!strcmp(argv[1], "B")) {
    type = BI;
  }
  else if (!strcmp(argv[1], "T")) {
    type = TAN;
  }
  else if (!strcmp(argv[1], "U")) {
    type = UNI;
  } else if (!strcmp(argv[1], "P")) {
    type = PSM;
  } else {
    throw invalid_argument("Choose a valid circuit type (U for univariate, B for bivariate and T for Tan et al.\n");
  }

  bool verbose = false;
  if (!strcmp(argv[8], "y"))
    verbose = true;

  //////////PARAMETER SET UP////////////////
  // Plaintext prime modulus
  unsigned long p = atol(argv[2]);
  // Field extension degree
  unsigned long d = atol(argv[3]);
  // Cyclotomic polynomial - defines phi(m)
  unsigned long m = atol(argv[4]);
  // Number of ciphertext prime bits in the modulus chain
  unsigned long nb_primes = atol(argv[5]);
  // Number of columns of Key-Switching matix (default = 2 or 3)
  unsigned long c = 3;

  if(type == PSM) {
    adjustingParameters(p, m, nb_primes, d);
    cout << "Parms: P " << p << " " << d << " " << m << " " << nb_primes << " " << argv[6] << " " << argv[7] << endl;
  }


  cout << "Initialising context object..." << endl;
  // Intialise context
  auto context = ContextBuilder<BGV>()
            .m(m)
            .p(p)
            .r(1)
            .bits(nb_primes)
            .c(c)
            .scale(6)
            .build();
  const EncryptedArray& ea = context.getEA();
  // Print the security level
  cout << "Ctx primes" << context.getCtxtPrimes() << endl;
  cout << "full primes" << context.fullPrimes() << endl;
  cout << "Q size: " << context.logOfProduct(context.getCtxtPrimes())/log(2.0) << endl;
  cout << "Q*P size: " << context.logOfProduct(context.fullPrimes())/log(2.0) << endl;
  cout << "Security: " << context.securityLevel() << endl;

  // Print the context
  context.getZMStar().printout();
  cout << endl;

  //maximal number of digits in a number
  unsigned long expansion_len = atol(argv[6]);

  // Secret key management
  // Create a secret key associated with the context
  SecKey secret_key(context);
  // Generate the secret key
  secret_key.GenSecKey();


  // Compute key-switching matrices that we need
  if(type == PSM) {

      const PAlgebra& al = ea.getPAlgebra();
      unsigned long slots = al.getNSlots();
      unsigned long enc_base = (p-1) >> 1;
      unsigned long maxsize = enc_base > slots ? slots : enc_base;

      for(uint g = 0; g < al.numOfGens(); g++  ) {
        for(uint r=1; r<maxsize; r <<= 1) {
          long v = al.coordinate(g, r);
          if( v!= 0) {
            secret_key.GenKeySWmatrix(1, context.getZMStar().genToPow(g, v), 0, 0);
          }
          v = al.coordinate(g, slots-r);
          if( v!= 0) {
            secret_key.GenKeySWmatrix(1, context.getZMStar().genToPow(g, v), 0, 0);
          }
        }
      }
      secret_key.setKeySwitchMap();

    //addSome1DMatrices(secret_key);

  } else if (expansion_len > 1) {
    if (context.getZMStar().numOfGens() == 1) {
      std::set<long> automVals;
      long e = 1;
      long ord = context.getZMStar().OrderOf(0);
      bool native = context.getZMStar().SameOrd(0);
      if(!native)
        automVals.insert(context.getZMStar().genToPow(0, -ord));
      while (e < expansion_len){
        long atm = context.getZMStar().genToPow(0, ord-e);
        automVals.insert(atm);
        e <<=1;
      }
      addTheseMatrices(secret_key, automVals);
    } else {
      addSome1DMatrices(secret_key);
    }
    if (d > 1 ) 
      addFrbMatrices(secret_key); //might be useful only when d > 1
  }

  // create Comparator (initialize after buildModChain)
  Comparator comparator(context, type, d, expansion_len, secret_key, verbose);

  //repeat experiments several times
  int runs = atoi(argv[7]);
  
  //test comparison circuit
  if(type == PSM) {
    comparator.test_compare_psm(runs);
  } else {
    comparator.test_compare(runs);
  }

  cout << " BS: " << static_cast<int>(log2((p - 1) >> 1)) << " S: " << context.securityLevel() << " - " << argv[0] << (type==PSM?" P ":" U ") << p << " " << d << " " << m << " " << nb_primes << " " << argv[6] << " " << argv[7] << endl;

  //printAllTimers(cout);

  return 0;
}
