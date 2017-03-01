// lfsr88.h: 32-bit Random number generator U[0,1): lfsr88
// Cycle length is (2^31 - 1)(2^29 - 1)(2^28 - 1) = 309485007947847626691444735 or approximately 2^88
// Author: Pierre L'Ecuyer,
// Source:  http://www.iro.umontreal.ca/~lecuyer/myftp/papers/tausme2.ps
// R. Saucier, 24 August 2016

#ifndef LFSR88_H
#define LFSR88_H

namespace LFSR88 {

static const uint32_t N_SEEDS = 3;
static const bitmatrix32_t MATRIX[N_SEEDS] = {
   {
      {
         0x00000000, 0x00002000, 0x00004000, 0x00008000, 0x00010000, 0x00020000, 0x00040001, 0x00080002,
         0x00100004, 0x00200008, 0x00400010, 0x00800020, 0x01000040, 0x02000080, 0x04000100, 0x08000200,
         0x10000400, 0x20000800, 0x40001000, 0x80000001, 0x00000002, 0x00000004, 0x00000008, 0x00000010,
         0x00000020, 0x00000040, 0x00000080, 0x00000100, 0x00000200, 0x00000400, 0x00000800, 0x00001000
      }
   },
   {
      {
         0x00000000, 0x00000000, 0x00000000, 0x00000080, 0x00000100, 0x00000200, 0x00000400, 0x00000800,
         0x00001000, 0x00002000, 0x00004000, 0x00008000, 0x00010000, 0x00020000, 0x00040000, 0x00080000,
         0x00100000, 0x00200000, 0x00400000, 0x00800000, 0x01000000, 0x02000000, 0x04000000, 0x08000001,
         0x10000002, 0x20000005, 0x4000000A, 0x80000014, 0x00000028, 0x00000050, 0x00000020, 0x00000040
      }
   },
   {
      {
         0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00200000, 0x00400000, 0x00800000, 0x01000000,
         0x02000001, 0x04000002, 0x08000004, 0x10000009, 0x20000012, 0x40000024, 0x80000048, 0x00000090,
         0x00000120, 0x00000240, 0x00000480, 0x00000900, 0x00001200, 0x00002400, 0x00004800, 0x00009000,
         0x00012000, 0x00024000, 0x00048000, 0x00090000, 0x00120000, 0x00040000, 0x00080000, 0x00100000
      }
   }
};
static const bitmatrix32_t MATRIX_INV[N_SEEDS] = {
   {
      {
         0x00000000, 0x00100000, 0x00200000, 0x00400000, 0x00800000, 0x01000000, 0x02000000, 0x04000000,
         0x08000000, 0x10000000, 0x20000000, 0x40000000, 0x80000001, 0x00000002, 0x00000004, 0x00000008,
         0x00000010, 0x00000020, 0x00000040, 0x00100080, 0x00200100, 0x00400200, 0x00800400, 0x01000800,
         0x02001000, 0x04002000, 0x08004000, 0x10008000, 0x20010000, 0x40020000, 0x80040000, 0x00080000
      }
   },
   {
      {
         0x00000000, 0x00000000, 0x00000000, 0x50000000, 0xa0000001, 0x40000002, 0x80000004, 0x00000008,
         0x00000010, 0x00000020, 0x00000040, 0x00000080, 0x00000100, 0x00000200, 0x00000400, 0x00000800,
         0x00001000, 0x00002000, 0x00004000, 0x00008000, 0x00010000, 0x00020000, 0x00040000, 0x00080000,
         0x00100000, 0x00200000, 0x00400000, 0x00800000, 0x01000000, 0x02000000, 0x54000000, 0xa8000000
      }
   },
   {
      {
         0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x49248000, 0x92490000, 0x24920000, 0x49240000,
         0x92480000, 0x24900000, 0x49200000, 0x92400000, 0x24800000, 0x49000000, 0x92000000, 0x24000000,
         0x48000000, 0x90000001, 0x20000002, 0x40000004, 0x80000008, 0x00000010, 0x00000020, 0x00000040,
         0x00000080, 0x00000100, 0x00000200, 0x00000400, 0x00000800, 0x49249000, 0x92492000, 0x24924000
      }
   }
};
static const uint32_t C0 = 0xfffffffful;   // 4294967295ul
static const uint32_t C1 = 0xfffffffeul;   // 4294967294ul
static const uint32_t C2 = 0xfffffff8ul;   // 4294967288ul
static const uint32_t C3 = 0xfffffff0ul;   // 4294967280ul

class lfsr88 : public Generator<uint32_t> {

public:
   lfsr88 ( void ) { // default constructor
   }

   lfsr88 ( std::vector<uint32_t> seed ) { // constructor from vector seed
   
      setState( seed );
   }
   
   virtual ~lfsr88() { // default destructor
   
      std::cout << "deleting lfsr88" << std::endl;
   }
   
   virtual void setState( std::vector<uint32_t> seed ) { // set the seeds
   
      assert( seed.size() >= N_SEEDS );
      
      // VERY IMPORTANT: The initial seeds _s[0], _s[1], _s[2] MUST be larger than 1, 7, and 15 respectively
      _s[0] = seed[0]; if ( _s[0] <  2 ) _s[0] +=  2;
      _s[1] = seed[1]; if ( _s[1] <  8 ) _s[1] +=  8;
      _s[2] = seed[2]; if ( _s[2] < 16 ) _s[2] += 16;
   }
   
   virtual void getState( std::vector<uint32_t>& seed ) { // get the seed vector
      
      assert( seed.size() >= N_SEEDS );
      for ( size_t i = 0; i < N_SEEDS; i++ ) seed[i] = _s[i];
   }
   
   virtual void jump_ahead( uintmax_t n ) {   // jump ahead the next n random numbers
      
      for ( size_t i = 0; i < N_SEEDS; i++ ) {
      
         Bitmatrix<uint32_t> A( MATRIX[i] );
         _s[i] = ( A^n ) * _s[i];
      }
   }
   
   virtual void jump_ahead( uintmax_t e, uintmax_t c ) {   // jump ahead the next n random numbers, where n = 2^e + c
   
      Bitmatrix<uint32_t> A, B;
      
      for ( size_t i = 0; i < N_SEEDS; i++ ) {
      
         if ( e ) {
            B = MATRIX[i];
            for ( uintmax_t j = 0; j < e; j++ ) B *= B;
         }
         A = MATRIX[i];
         A = A^c;
         if ( e ) A *= B;
         _s[i] = A * _s[i];
      }
   }
   
   virtual void jump_back( uintmax_t n ) {   // jump ahead the next n random numbers
      
      for ( size_t i = 0; i < N_SEEDS; i++ ) {
      
         Bitmatrix<uint32_t> A( MATRIX_INV[i] );
         _s[i] = ( A^n ) * _s[i];
      }
   }
   
   virtual void jump_back( uintmax_t e, uintmax_t c ) {   // jump ahead the next n random numbers, where n = 2^e + c
   
      Bitmatrix<uint32_t> A, B;
      
      for ( size_t i = 0; i < N_SEEDS; i++ ) {
      
         if ( e ) {
            B = MATRIX_INV[i];
            for ( uintmax_t j = 0; j < e; j++ ) B *= B;
         }
         A = MATRIX_INV[i];
         A = A^c;
         if ( e ) A *= B;
         _s[i] = A * _s[i];
      }
   }
   
   virtual void jump_cycle( void ) {
   
      const uint32_t A = 31, B = 29, C = 28;
      jump_ahead( A + B + C, 0 );
      jump_back( A + B, 0 );
      jump_back( A + C, 0 );
      jump_back( B + C, 0 );
      jump_ahead( A, 0 );
      jump_ahead( B, 0 );
      jump_ahead( C, 0 );
      jump_back( 1 );
   }
/*
	virtual void jump_cycle( void ) { // jump ahead a full cycle of lfsr88
	
	   std::bitset<88> p( std::string( "1111111111111111111111111110011000000000000000000000000010101111111111111111111111111111" ) );
      for ( size_t i = 0; i < p.size(); ++i ) if ( p.test(i) ) jump_ahead( i, 0 );
	}
*/
   
   virtual uint32_t rng32( void ) {   // returns 32-bit integer
   
      _s[0] = ( ( _s[0] & C1 ) << 12 ) ^ ( ( ( _s[0] << 13 ) ^ _s[0] ) >> 19 ) ;
      _s[1] = ( ( _s[1] & C2 ) <<  4 ) ^ ( ( ( _s[1] <<  2 ) ^ _s[1] ) >> 25 ) ;
      _s[2] = ( ( _s[2] & C3 ) << 17 ) ^ ( ( ( _s[2] <<  3 ) ^ _s[2] ) >> 11 ) ;
   
      return ( _s[0] ^ _s[1] ^ _s[2] ) & C0;
   }

   virtual uint64_t rng64( void ) {   // returns 64-bit integer
   
      uint64_t low  = rng32();
      uint64_t high = rng32();
      return low | ( high << 32 );
   }

   virtual double rng32_01( void ) {   // returns a double in [0,1)
   
      return double( rng32() ) * TWO32_INV;
   }

   virtual long double rng64_01( void ) {   // returns a long double in [0,1)
   
      return double( rng64() ) * TWO64_INV;
   }

private:
   
   uint32_t _s[ N_SEEDS ];
   
}; // end lfsr88 class
} // end namespace LFSR88

#endif // LFSR88_H
