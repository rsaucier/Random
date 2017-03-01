// jkiss.h: Based upon Marsaglia's Keep It Simple Stupid RNG
// Ref: Good Practice in (Pseudo) Random Number Generation for Bioinformatics Applications
//      David Jones, UCL Bioinformatics Group (d.jones@cs.ucl.ac.uk), May 7, 2010
// Period is (2^32)(2^32-1)(4294584393(2^31)-1) = 17012601507030308243410262827431100416 or approximately 2^127
// Period of MWC is 4294584393(2^31)-1 = 9222549758923505663, where 4294584393 = 2^32 - 2^19 + 2^17 + 2^13 + 2^11 + 2^6 + 2^3 + 1
// R. Saucier, 24 August 2016

#ifndef JKISS_H
#define JKISS_H

namespace JKISS {
   
   static const bitmatrix32_t MATRIX = { // 32x32 bitmatrix
      {
         0x08400021, 0x10800042, 0x21400085, 0x4280010a, 0x85000214, 0x0a000428, 0x14000850, 0x284010a1, 
         0x50802142, 0xa1004284, 0x42008508, 0x84010a10, 0x08021420, 0x10042840, 0x20085080, 0x4010a100, 
         0x80214200, 0x00428400, 0x00850800, 0x010a1000, 0x02142000, 0x04284000, 0x08508000, 0x10a10000, 
         0x21420000, 0x42840000, 0x85080000, 0x08100000, 0x10200000, 0x20400000, 0x40800000, 0x81000000
      }
   };
   static const bitmatrix32_t MATRIX_INV = { // A_INV = A + A^7 + A^9 + A^10 + A^11 + A^13 + A^19 + A^20 + A^21 + A^22 + A^23 + A^31, where A = MATRIX
      {
         0x9ce52d63, 0x39ca5ac6, 0x7394b58c, 0xe7296b18, 0xce52d630, 0x9ca5ac60, 0x7b5bdce1, 0xb4a73de3,
         0x694e7bc6, 0xd29cf78c, 0x5294a508, 0xa5294a10, 0x4a529420, 0x94a52840, 0x6b5ad4a1, 0xd6b5a942,
         0xad6b5284, 0x5ad6a508, 0xb5ad4a10, 0x6b5a9420, 0xd6b52840, 0xef7ad4a1, 0xdef5a942, 0xbdeb5284,
         0x7bd6a508, 0xf7ad4a10, 0xef5a9420, 0xdeb52840, 0xff7ad4a1, 0xfef5a942, 0xfdeb5284, 0xfbd6a508
      }
   };
   static const uint32_t LC_MULT      = 0x12bf507dul;            // 314527869ul;
   static const uint32_t LC_CONST     = 0x0012d687ul;            // 1234567ul;
   static const uint32_t LC_MULT_INV  = 0x6200a8d5ul;            // 1644210389ul;
   static const uint64_t MWC_MULT     = 0x00000000fffa2849ull;   // 4294584393ull;
   static const uint64_t MWC_MOD      = 0xfffa2848ffffffffull;   // 18445099517847011327ull;
   static const uint64_t MWC_MULT_INV = 0x0000000100000000ull;   // 4294967296ull;
   static const uint64_t LC_PERIOD    = 0x0000000100000000ull;   // 4294967296ull;
   static const uint64_t SR_PERIOD    = 0x00000000ffffffffull;   // 4294967295ull;
   static const uint64_t MWC_PERIOD   = 0x7ffd14247fffffffull;   // 9222549758923505663ull;
   static const uint32_t N_SEEDS      = 4;                       // requires four 32-bit words

class jkiss : public Generator<uint32_t> {

public:
   jkiss( void ) { // default constructor
   }
   
   jkiss( std::vector<uint32_t> seed ) { // constructor from seed vector
   
      setState( seed );
   }
   
   virtual ~jkiss() { // default destructor
   
      std::cout << "deleting jkiss" << std::endl;
   }
   
   virtual void setState( std::vector<uint32_t> seed ) { // set the seeds
   
      assert( seed.size() >= N_SEEDS );
      _s1 = seed[0];
      _s2 = seed[1];
      _s3 = seed[2];
      _s4 = seed[3];
   }
   
   virtual void getState( std::vector<uint32_t>& seed ) { // get the seed vector
   
      assert( seed.size() >= N_SEEDS );
      seed[0] = _s1;
      seed[1] = _s2;
      seed[2] = _s3;
      seed[3] = _s4;
   }
   
   virtual void jump_ahead( uintmax_t n ) { // jump ahead the next n random numbers

#ifdef UINT64_C // use the native 64-bit capability
   
      uint64_t p = mul_mod64( pow_mod64( LC_MULT, n, M ), _s1, M );
      uint64_t q = mul_mod64( LC_CONST, gs_mod64( LC_MULT, n, M ), M );
      _s1 = static_cast<uint32_t>( add_mod64( p, q, M ) );
   
#else // native 64-bit not available, so use double instead
   
      double p = mul_mod( pow_mod( LC_MULT, n, M ), _s1, M );
      double q = mul_mod( LC_CONST, gs_mod( LC_MULT, n, M ), M );
      _s1 = static_cast<uint32_t>( add_mod( p, q, M ) );
   
#endif // UINT64_C

      Bitmatrix<uint32_t> A, B( MATRIX );
      A = B^n;
      _s2 = A * _s2;
      
      uint64_t a = _s3 + ( (uint64_t)_s4 << 32u );
      a = mul_mod64( pow_mod64( MWC_MULT, n, MWC_MOD ), a, MWC_MOD );
      _s4 = ( uint32_t )( a >> 32u );
      _s3 = ( uint32_t )( a );
   }
   
   virtual void jump_ahead( uintmax_t e, uintmax_t c ) {   // jump ahead the next n random numbers, where n = 2^e + c
      
#ifdef UINT64_C // use the native 64-bit capability

      uint64_t p = mul_mod64( pow_mod64( LC_MULT, e, c, M ), _s1, M );
      uint64_t q = mul_mod64( LC_CONST, gs_mod64( LC_MULT, e, c, M ), M );
      _s1 = static_cast<uint32_t>( add_mod64( p, q, M ) );

#else // native 64-bit not available, so use double instead
   
      double p = mul_mod( pow_mod( LC_MULT, e, c, M ), _s1, M );
      double q = mul_mod( LC_CONST, gs_mod( LC_MULT, e, c, M ), M );
      _s1 = static_cast<uint32_t>( add_mod( p, q, M ) );
   
#endif // UINT64_C

      Bitmatrix<uint32_t> A, B;
   
      if ( e ) {
         B = MATRIX;
         for ( size_t i = 0; i < e; i++ ) B *= B;
      }
      A = MATRIX;
      A = A^c;
      if ( e ) A *= B;
      _s2 = A * _s2;
      
      uint64_t a = _s3 + ( (uint64_t)_s4 << 32u );
      a = mul_mod64( pow_mod64( MWC_MULT, e, c, MWC_MOD ), a, MWC_MOD );
      _s4 = ( uint32_t )( a >> 32u );
      _s3 = ( uint32_t )( a );
   }
   
   virtual void jump_back( uintmax_t n ) { // jump back the next n random numbers
   
#ifdef UINT64_C // use the native 64-bit capability
   
      uint64_t p = mul_mod64( pow_mod64( LC_MULT_INV, n, M ), add_mod64( _s1, -LC_CONST, M ), M );
      uint64_t q = mul_mod64( -LC_CONST, gs_mod64( LC_MULT_INV, n, M ), M );
      uint64_t r = add_mod64( p, q, M );
      _s1 = static_cast<uint32_t>( add_mod( LC_CONST, r, M ) );

#else // native 64-bit not available, so use double instead
   
      double p = mul_mod( pow_mod( LC_MULT_INV, n, M ), add_mod( _s1, -LC_CONST, M ), M );
      double q = mul_mod( -LC_CONST, gs_mod( LC_MULT_INV, n, M ), M );
      double r = add_mod( p, q, M );
      _s1 = static_cast<uint32_t>( add_mod( LC_CONST, r, M ) );
   
#endif // UINT64_C

      Bitmatrix<uint32_t> A( MATRIX_INV );
      A = A^n;
      _s2 = A * _s2;

      uint64_t a = _s3 + ( (uint64_t)_s4 << 32u );
      a = mul_mod64( pow_mod64( MWC_MULT_INV, n, MWC_MOD ), a, MWC_MOD );
      _s4 = ( uint32_t )( a >> 32u );
      _s3 = ( uint32_t )( a );
   }
   
   virtual void jump_back( uintmax_t e, uintmax_t c ) {   // jump back the next n random numbers, where n = 2^e + c
   
#ifdef UINT64_C // use the native 64-bit capability
   
      uint64_t p = mul_mod64( pow_mod64( LC_MULT_INV, e, c, M ), add_mod64( _s1, -LC_CONST, M ), M );
      uint64_t q = mul_mod64( -LC_CONST, gs_mod64( LC_MULT_INV, e, c, M ), M );
      uint64_t r = add_mod64( p, q, M );
      _s1 = static_cast<uint32_t>( add_mod64( LC_CONST, r, M ) );
   
#else // native 64-bit not available, so use double instead
   
      double p = mul_mod( pow_mod( LC_MULT_INV, e, c, M ), add_mod( _s1, -LC_CONST, M ), M );
      double q = mul_mod( -LC_CONST, gs_mod( LC_MULT_INV, e, c, M ), M );
      double r = add_mod( p, q, M );
      _s1 = static_cast<uint32_t>( add_mod( LC_CONST, r, M ) );
   
#endif // UINT64_C

      Bitmatrix<uint32_t> A, B;
   
      if ( e ) {
         B = MATRIX_INV;
         for ( size_t i = 0; i < e; i++ ) B *= B;
      }
      A = MATRIX_INV;
      A = A^c;
      if ( e ) A *= B;
      _s2 = A * _s2;
      
      uint64_t a = _s3 + ( (uint64_t)_s4 << 32u );
      a = mul_mod64( pow_mod64( MWC_MULT_INV, e, c, MWC_MOD ), a, MWC_MOD );
      _s4 = ( uint32_t )( a >> 32u );
      _s3 = ( uint32_t )( a );
   }
   
   virtual void jump_cycle( void ) { // jump ahead a full cycle of jkiss

      jump_ahead( 127, 0 );
      jump_back( 114, 0 );
      jump_ahead( 112, 0 );
      jump_ahead( 108, 0 );
      jump_ahead( 106, 0 );
      jump_ahead( 101, 0 );
      jump_ahead( 98, 0 );
      jump_ahead( 82, 0 );
      jump_back( 80, 0 );
      jump_back( 76, 0 );
      jump_back( 74, 0 );
      jump_back( 69, 0 );
      jump_back( 66, 0 );
      jump_back( 64, 0 );
      jump_back( 63, 0 );
      jump_ahead( 32, 0 );
   }
/*
	virtual void jump_cycle( void ) { // jump ahead a full cycle of jkiss
	
	   std::bitset<127> p( std::string( "1111111111111010001010000100100000000000000001011101011110110101000000000000000000000000000000100000000000000000000000000000000" ) );
      for ( size_t i = 0; i < p.size(); ++i ) if ( p.test(i) ) jump_ahead( i, 0 );
	}
*/
   virtual uint32_t rng32( void ) { // returns the next random number (as a 32-bit unsigned int)
   
      _s1 = LC_MULT * _s1 + LC_CONST;
   
      _s2 ^= ( _s2 << 5 ), _s2 ^= ( _s2 >> 7 ), _s2 ^= ( _s2 << 22 );
   
      uint64_t a = MWC_MULT * _s3 + _s4;
      _s4 = ( a >> 32u );
      _s3 = uint32_t( a );
   
      return _s1 + _s2 + _s3;
   }
   
   virtual uint64_t rng64( void ) {   // returns 64-bit unsigned integer
   
      uint64_t low  = rng32();
      uint64_t high = rng32();
      return low | ( high << 32 );
   }
   
   virtual double rng32_01( void ) { // returns a double in the half-open interval [0,1)
   
      return double( rng32() ) * TWO32_INV;
   }
   
   virtual long double rng64_01( void ) {   // returns a long double in [0,1)
   
      return double( rng64() ) * TWO64_INV;
   }
   
private:

   uint32_t _s1, _s2, _s3, _s4;

}; // end jkiss class
} // end namespace JKISS

#endif // JKISS_H
