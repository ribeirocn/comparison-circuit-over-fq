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
* This code implements the Private Set Membership primitive described in [Ribeiro23].
* The code was adapted from https://github.com/iliailia/comparison-circuit-over-fq
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


// the main function that takes 7 arguments (type in Terminal: ./psm_circuit argv[1] argv[2] argv[3] argv[4] argv[5] argv[6] argv[7] argv[8])
// argv[1] - circuit type (I - Integer or S - String)
// argv[2] - the plaintext modulus
// argv[3] - the dimension of a vector space over a finite field
// argv[4] - the order of the cyclotomic ring
// argv[5] - the bitsize of the ciphertext modulus in ciphertexts (HElib increases it to fit the moduli chain). The modulus used for public-key generation
// argv[6] - the length of vectors to be compared
// argv[7] - the number of strings to be compared
// argv[8] - the number of experiment repetitions
// argv[9] - print debug info (y/n)

// some parameters for quick testing
// String comparasion with UniSlot packing
// S 257 16 31523 480 1 1000 1 y
// S 521 16 37193 580 1 1000 1 y
// S 1031 16 32743 580 1 1500 1 y
// S 65537 16 74789 950 1 500 1 y

// String comparasion with MultiSlot packing
// S 257 1 31523 480 16 90 1 y
// S 521 1 36517 580 16 100 1 Y
// S 1031 1 32743 580 16 100 1 y
// S 65537 1 74703 950 16 500 1 y

void adjustingParameters(unsigned long& p, unsigned long& m, unsigned long nb_primes, unsigned long d, long ss_size=-1);

int main(int argc, char *argv[])
{
  if (argc < 10)
  {
    throw invalid_argument("There should be exactly 9 arguments\n");
  }

  CircuitType type = UNI;
  if (!strcmp(argv[1], "I"))
  {
    type = PSM;
  }
  else if (!strcmp(argv[1], "S"))
  {
    type = PSMS;
  }
  else
  {
    throw invalid_argument("Choose a valid circuit type (S for String, I for Integer\n");
  }

  bool verbose = false;
  if (!strcmp(argv[9], "y"))
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

  // maximal number of digits in a number
  unsigned long expansion_len = atol(argv[6]);
  unsigned long ss_size = atol(argv[7]);

  adjustingParameters(p, m, nb_primes, d, (long)expansion_len*ss_size);
  cout << "Parms: S " << p << " " << d << " " << m << " " << nb_primes << " " << argv[6] << " " << argv[7] << " " << argv[8] << endl;

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
  const EncryptedArray &ea = context.getEA();
  // Print the security level
  cout << "Ctx primes" << context.getCtxtPrimes() << endl;
  cout << "full primes" << context.fullPrimes() << endl;
  cout << "Q size: " << context.logOfProduct(context.getCtxtPrimes()) / log(2.0) << endl;
  cout << "Q*P size: " << context.logOfProduct(context.fullPrimes()) / log(2.0) << endl;
  cout << "Security: " << context.securityLevel() << endl;

  // Print the context
  context.getZMStar().printout();
  cout << endl;



  // Secret key management
  // Create a secret key associated with the context
  SecKey secret_key(context);
  // Generate the secret key
  secret_key.GenSecKey();

  const PAlgebra &al = ea.getPAlgebra();
  unsigned long slots = al.getNSlots();
  unsigned long enc_base = (p - 1) >> 1;
  unsigned long maxsize = enc_base > slots ? slots : enc_base;

  for (uint g = 0; g < al.numOfGens(); g++) {
    for (uint r = 1; r < maxsize; r <<= 1) {
      long v = al.coordinate(g, r);
      if (v != 0) {
        secret_key.GenKeySWmatrix(1, context.getZMStar().genToPow(g, v), 0, 0);
      }
      v = al.coordinate(g, slots - r);
      if (v != 0) {
        secret_key.GenKeySWmatrix(1, context.getZMStar().genToPow(g, v), 0, 0);
      }
    }
  }
  secret_key.setKeySwitchMap();

  // addSome1DMatrices(secret_key);

  if (d > 1 )
    addFrbMatrices(secret_key); // might be useful only when d > 1

  // create Comparator (initialize after buildModChain)
  Comparator comparator(context, type, d, expansion_len, secret_key, verbose, ss_size);

  // repeat experiments several times
  int runs = atoi(argv[8]);

  // test comparison circuit
  comparator.test_string_psm(runs);

  cout << " SS: " << argv[7] << " S: " << context.securityLevel() << " - " << argv[0] << " " << argv[1] << " " << p << " " << d << " " << m << " " << nb_primes << " " << argv[6] << " " << argv[7] << " " << argv[8] << endl;

  // printAllTimers(cout);

  return 0;
}
