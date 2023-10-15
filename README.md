# Faster homomorphic comparison operations for BGV and BFV extended for comparasion with the PSM

This repository contains the C++ code of the homomorphic pattern matching algorithm from the paper ["Faster homomorphic comparison operations for BGV and BFV"](https://eprint.iacr.org/2021/315) by Ilia Iliashenko and Vincent Zucca, extended with the PSM primitive for comparison 

## Installation guide
To build the code, install [HElib](https://github.com/homenc/HElib) (>= 2.0.0), 

    mkdir build
    cd build
    cmake ../src
    make

The binaries will be available at build/bin

## How to use
### Integer comparison
To test the basic comparison of integers, use the following command
  
    ./comparison_circuit circuit_type p d m q l runs print_debug_info
    
where
+ `circuit_type` takes one of four values `P`, `U`, `B` or `T`. Th first one corresponds to our Private Set Membership Primitive, and the remaining three corresponde to the univariate, bivariate circuits from ["Faster homomorphic comparison operations for BGV and BFV"](https://eprint.iacr.org/2021/315), and the circuit of [Tan et al.](https://eprint.iacr.org/2019/332).
+ `p`: the plaintext modulus, must be a prime number.
+ `d`: the dimension of a vector space over the slot finite field.
+ `m`: the order of the cyclotomic ring.
+ `q`: the minimal bitsize of the ciphertext modulus in ciphertexts. The actual size of the modulus is automatically chosen by HElib.
+ `l`: the length of finite field vectors to be compared.
+ `runs`: the number of experiments.
+ `print_debug_info`: type `y` or `n` to show/hide more details on computation.
More details on these parameters can be found in Section 5 of the paper.

The following lines compares the execution of a PSM test with a univariate 1, when comparing two numbers of 15 bits over a ring of order 65336, with length l. Notice that we are only comparing two numbers which is not the best scenario for the univariate scheme.
  
    ./comparison_circuit P 65537 1 65536 730 1 10 y
    ./comparison_circuit U 65537 1 65536 730 1 1 y
### String Comparison
String comparison is different in the sense that two strings are naturally divided in digits of character size
Strings may be packing in two different ways. One more compact and one more efficient. The UniSlot is more compact and the MultiSlot is more efficient

    ./psm_circuit S p d m q l N runs print_debug_info
    
where the most arguments are analogous to the comparison circuit above. The additional argument is
+ `N`: the number of strings comprising the set to search.

The size of the string to compared is defined by d*l. The UniSlot vs MultiSlot is defined by parameters l and d. When l=1 then d>1 and is the UniSlot packing, when d=1 the l>1 and it is the Multislot packing.

Two running exemples are
    ./psm_circuit S 257 1 31523 480 16 90 1 y
    ./psm_circuit S 257 16 31523 480 1 1000 1 y