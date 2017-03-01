// Generator.h: template class file for random number generators
// R. Saucier, July 2016

#ifndef GENERATOR_H
#define GENERATOR_H

#include "Bitmatrix.h"
#include "mod_math.h"
#include <vector>
#include <bitset>
#include <iostream>

template <class T>   // for 32-bit and 64-bit generators
class Generator {

public:

   virtual ~Generator() {};// std::cout << "deleting Generator" << std::endl; }
   virtual void setState( std::vector<T> seed ) = 0;
   virtual void getState( std::vector<T>& seed ) = 0;
   virtual void jump_ahead( uintmax_t ) = 0;
   virtual void jump_ahead( uintmax_t,  uintmax_t ) = 0;
   virtual void jump_back( uintmax_t ) = 0;
   virtual void jump_back( uintmax_t,  uintmax_t ) = 0;
   virtual void jump_cycle( void ) = 0;

   virtual uint32_t    rng32( void ) = 0;      // returns 32-bit integer
   virtual uint64_t    rng64( void ) = 0;      // returns 64-bit integer
   virtual double      rng32_01( void ) = 0;   // returns double in [0,1)
   virtual long double rng64_01( void ) = 0;   // returns long double in [0,1)

   inline double u32( double a = 0., double b = 1. ) { return a + ( b - a ) * this->rng32_01(); }
   inline double u64( double a = 0., double b = 1. ) { return a + ( b - a ) * this->rng64_01(); }
};

// 32-bit generators
#include "kiss.h"
#include "jkiss.h"
#include "lfsr88.h"
#include "lfsr113.h"

// 64-bit generators
#include "jlkiss.h"
#include "jlkiss64.h"
#include "lfsr258.h"

#endif
