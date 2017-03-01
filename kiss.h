// kiss.h: Marsaglia's Keep It Simple Stupid RNG, which consists of a combination of linear congruential, 3-shift-register, and multiply with carry
// Period is (2^32)(2^32-1)(698769069(2^31)-1) = 27681094672891588090390813844460011520, approximately 2^124,
// where 698769069 = 2^29 + 2^27 + 2^24 + 2^23 + 2^21 + 2^18 + 2^17 + 2^14 + 2^12 + 2^11 + 2^10 + 2^9 + 2^7 + 2^5 + 2^3 + 2^2 + 1
// R. Saucier, 24 August 2016

#ifndef KISS_H
#define KISS_H
#include <bitset>

namespace KISS {

   static const bitmatrix32_t MATRIX = {
      {
         0x00042021, 0x00084042, 0x00108084, 0x00210108, 0x00420231, 0x00840462, 0x010808C4, 0x02101188,
         0x04202310, 0x08404620, 0x10808C40, 0x21011880, 0x42023100, 0x84046200, 0x0808C400, 0x10118800,
         0x20231000, 0x40462021, 0x808C4042, 0x01080084, 0x02100108, 0x04200210, 0x08400420, 0x10800840,
         0x21001080, 0x42002100, 0x84004200, 0x08008400, 0x10010800, 0x20021000, 0x40042000, 0x80084000
      }
   };
   static const bitmatrix32_t MATRIX_INV = { 
      {
         0xf2b58529, 0xe56b0a52, 0xded6b4a5, 0xbdad694a, 0x7b5ad294, 0xf6b5a528, 0xed6b4a50, 0xced634a1,
         0x9dac6942, 0x3b58d284, 0x76b1a508, 0xed634a10, 0xcec63421, 0x9d8c6842, 0x3b18d084, 0x7631a108,
         0xec634210, 0xccc62421, 0x998c4842, 0x33189084, 0x66312108, 0xcc624210, 0x88c40420, 0x11880840,
         0x23101080, 0x46202100, 0x8c404200, 0x08800400, 0x11000800, 0x22001000, 0x44002000, 0x88004000
      }
   };
   static const uint32_t LC_MULT         = 0x00010dcd;           // 69069UL
   static const uint32_t LC_CONST        = 0x00003039;           // 12345UL
   static const uint32_t LC_MULT_INV     = 0xa5e2a705;           // 2783094533UL
   static const uint64_t MWC_MULT        = 0x0000000029a65ead;   // 698769069ull;
   static const uint64_t MWC_MOD         = 0x29a65eacffffffff;   // 3001190298811367423ULL
   static const uint64_t MWC_MULT_INV    = 0x0000000100000000;   // 4294967296ULL;
   static const uint32_t N_SEEDS         = 4;                    // requires four 32-bit seeds

class kiss : public Generator<uint32_t> {

public:
   kiss( void ) { // default constructor
   }
   
   kiss( std::vector<uint32_t> seed ) { // constructor from seed vector
   
      setState( seed );
   }
   
   virtual ~kiss() {   // default destructor
   
      std::cout << "deleting kiss" << std::endl;
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
   
   virtual void jump_ahead( uintmax_t n ) { // jumps ahead the next n random numbers
      
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
   
   virtual void jump_back( uintmax_t n ) { // jump back n 
      
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
/*
   virtual void jump_cycle( void ) { // jump ahead a full cycle of kiss

      jump_ahead( 124, 0 );
      jump_ahead( 122, 0 );
      jump_ahead( 120, 0 );
      jump_back( 118, 0 );
      jump_ahead( 116, 0 );
      jump_ahead( 114, 0 );
      jump_back( 112, 0 );
      jump_ahead( 110, 0 );
      jump_back( 109, 0 );
      jump_ahead( 108, 0 );
      jump_back( 104, 0 );
      jump_ahead( 102, 0 );
      jump_ahead( 100, 0 );
      jump_ahead( 99, 0 );
      jump_back( 97, 0 );
      jump_ahead( 95, 0 );
      jump_back( 93, 0 );
      jump_ahead( 91, 0 );
      jump_ahead( 90, 0 );
      jump_back( 88, 0 );
      jump_ahead( 85, 0 );
      jump_ahead( 84, 0 );
      jump_back( 82, 0 );
      jump_ahead( 80, 0 );
      jump_back( 78, 0 );
      jump_ahead( 76, 0 );
      jump_ahead( 71, 0 );
      jump_ahead( 69, 0 );
      jump_ahead( 67, 0 );
      jump_ahead( 63, 0 );
      jump_ahead( 32, 0 );
   }
*/
	virtual void jump_cycle( void ) { // jump ahead a full cycle of kiss
	
	   std::bitset<124> p( std::string( "10100110100110010111101010110011010110010110011010000101010001000000000000000000000000000000100000000000000000000000000000000" ) );
      for ( size_t i = 0; i < p.size(); ++i ) if ( p.test(i) ) jump_ahead( i, 0 );
	}

   virtual uint32_t rng32( void ) { // returns the next random number (as a 32-bit unsigned int)
   
      _s1 = LC_MULT * _s1 + LC_CONST;
   
      _s2 ^= ( _s2 << 13 ), _s2 ^= ( _s2 >> 17 ), _s2 ^= ( _s2 << 5 );
   
      uint64_t a = MWC_MULT * _s3 + _s4;
      _s4 = ( a >> 32u );
      _s3 = uint32_t( a );
   
      return _s1 + _s2 + _s3;
   }
   
   virtual uint64_t rng64( void ) {   // returns 64-bit integer
   
      uint64_t low  = rng32();
      uint64_t high = rng32();
      return low | ( high << 32 );
   }
   
   virtual double rng32_01( void ) { // returns a random number in the half-open interval [0,1)
   
      return double( rng32() ) * TWO32_INV;
   }
   
   virtual long double rng64_01( void ) {   // returns a long double in [0,1)
   
      return double( rng64() ) * TWO64_INV;
   }
   
private:

   uint32_t _s1, _s2, _s3, _s4;

}; // end kiss class
} // end namespace KISS

#endif // KISS_H
