// 32-bit Random number generator U[0,1): lfsr113
// Period is (2^31 - 1)(2^29 - 1)(2^28 - 1)(2^25 - 1) = 10384593344720504788331840650870785 or approximately 2^113
// Author: Pierre L'Ecuyer,
// Source: http://www.iro.umontreal.ca/~lecuyer/myftp/papers/tausme2.ps
// R. Saucier, 24 August 2016

#ifndef LFSR113_H
#define LFSR113_H

namespace LFSR113 {

static const uint32_t N_SEEDS = 4;
static const bitmatrix32_t MATRIX[N_SEEDS] = {
   {
      {
         0x00000000, 0x00080000, 0x00100000, 0x00200000, 0x00400000, 0x00800000, 0x01000000, 0x02000001,
         0x04000002, 0x08000004, 0x10000008, 0x20000010, 0x40000020, 0x80000041, 0x00000082, 0x00000104,
         0x00000208, 0x00000410, 0x00000820, 0x00001040, 0x00002080, 0x00004100, 0x00008200, 0x00010400,
         0x00020800, 0x00041000, 0x00002000, 0x00004000, 0x00008000, 0x00010000, 0x00020000, 0x00040000
      }
   },
	{
      {
         0x00000000, 0x00000000, 0x00000000, 0x00000020, 0x00000040, 0x00000080, 0x00000100, 0x00000200,
         0x00000400, 0x00000800, 0x00001000, 0x00002000, 0x00004000, 0x00008000, 0x00010000, 0x00020000,
         0x00040000, 0x00080000, 0x00100000, 0x00200000, 0x00400000, 0x00800000, 0x01000000, 0x02000000,
         0x04000000, 0x08000001, 0x10000002, 0x20000005, 0x4000000a, 0x80000014, 0x00000008, 0x00000010
      }
   },
	{
      {
         0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000800, 0x00001000, 0x00002000, 0x00004000,
         0x00008001, 0x00010002, 0x00020004, 0x00040008, 0x00080010, 0x00100020, 0x00200040, 0x00400080,
         0x00800100, 0x01000200, 0x02000400, 0x04000000, 0x08000000, 0x10000001, 0x20000002, 0x40000004,
         0x80000008, 0x00000010, 0x00000020, 0x00000040, 0x00000080, 0x00000100, 0x00000200, 0x00000400
      }
   },
	{
      {
         0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00100000,
         0x00200000, 0x00400001, 0x00800002, 0x01000004, 0x02000009, 0x04000012, 0x08000024, 0x10000048,
         0x20000090, 0x40000120, 0x80000240, 0x00000480, 0x00000900, 0x00001200, 0x00002400, 0x00004800,
         0x00009000, 0x00012000, 0x00024000, 0x00048000, 0x00090000, 0x00020000, 0x00040000, 0x00080000
      }
   }
};
static const bitmatrix32_t MATRIX_INV[N_SEEDS] = {
   {
      {
         0x00000000, 0x04104000, 0x08208000, 0x10410000, 0x20820000, 0x41040000, 0x82080000, 0x04100000,
         0x08200000, 0x10400000, 0x20800000, 0x41000000, 0x82000000, 0x04000000, 0x08000000, 0x10000000,
         0x20000000, 0x40000000, 0x80000001, 0x00000002, 0x00000004, 0x00000008, 0x00000010, 0x00000020,
         0x00000040, 0x00000080, 0x04104100, 0x08208200, 0x10410400, 0x20820800, 0x41041000, 0x82082000
      }
   },
   {
      {
         0x00000000, 0x00000000, 0x00000000, 0x40000002, 0x80000004, 0x00000008, 0x00000010, 0x00000020,
         0x00000040, 0x00000080, 0x00000100, 0x00000200, 0x00000400, 0x00000800, 0x00001000, 0x00002000,
         0x00004000, 0x00008000, 0x00010000, 0x00020000, 0x00040000, 0x00080000, 0x00100000, 0x00200000,
         0x00400000, 0x00800000, 0x01000000, 0x02000000, 0x04000000, 0x08000001, 0x50000000, 0xa0000001
      }
   },
   {
      {
         0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x02000000, 0x04000000, 0x08000000, 0x10000001,
         0x20000002, 0x40000004, 0x80000008, 0x00000010, 0x00000020, 0x00000040, 0x00000080, 0x00000100,
         0x00000200, 0x00000400, 0x00000800, 0x02001000, 0x04002000, 0x08004000, 0x10008000, 0x20010000,
         0x40020000, 0x80040000, 0x00080000, 0x00100000, 0x00200000, 0x00400000, 0x00800000, 0x01000000
      }
   },
   {
      {
         0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x92480000,
         0x24900000, 0x49200000, 0x92400000, 0x24800000, 0x49000000, 0x92000001, 0x24000002, 0x48000004,
         0x90000008, 0x20000010, 0x40000020, 0x80000040, 0x00000080, 0x00000100, 0x00000200, 0x00000400,
         0x00000800, 0x00001000, 0x00002000, 0x00004000, 0x00008000, 0x92490000, 0x24920000, 0x49240000
      }
   }
};
static const uint32_t C0 = 0xfffffffful;   // 4294967295ul
static const uint32_t C1 = 0xfffffffeul;   // 4294967294ul
static const uint32_t C2 = 0xfffffff8ul;   // 4294967288ul
static const uint32_t C3 = 0xfffffff0ul;   // 4294967280ul
static const uint32_t C4 = 0xffffff80ul;   // 4294967168ul

class lfsr113 : public Generator<uint32_t> {

public:
   lfsr113 ( void ) { // default constructor
   }

   lfsr113 ( std::vector<uint32_t> seed ) { // constructor from vector seed
   
      setState( seed );
   }
   
   virtual ~lfsr113() { // default destructor
	
      std::cout << "deleting lfsr113" << std::endl;
   }
   
   virtual void setState( std::vector<uint32_t> seed ) { // set the seeds
   
      assert( seed.size() >= N_SEEDS );
		
      // VERY IMPORTANT: The initial seeds _s1, _s2, _s3 MUST be larger than 1, 7, 15, and 127 respectively
      _s[0] = seed[0]; if ( _s[0] <   2 ) _s[0] +=   2;
      _s[1] = seed[1]; if ( _s[1] <   8 ) _s[1] +=   8;
      _s[2] = seed[2]; if ( _s[2] <  16 ) _s[2] +=  16;
      _s[3] = seed[3]; if ( _s[3] < 128 ) _s[3] += 128;
   }
   
   virtual void getState( std::vector<uint32_t>& seed ) { // get the seed vector
      
		assert( seed.size() >= N_SEEDS );
      for ( size_t i = 0; i < N_SEEDS; i++ ) seed[i] = _s[i];
   }
   
   virtual void jump_ahead( uintmax_t n ) { // jumps ahead the next n random numbers
		
		for ( size_t i = 0; i < N_SEEDS; i++ ) {
      
         Bitmatrix<uint32_t> A( MATRIX[i] );
         _s[i] = ( A^n ) * _s[i];
      }
   }
   
   virtual void jump_ahead( uintmax_t e, uintmax_t c ) {   // jumps ahead the next n random numbers, where n = 2^e + c
		
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
   
   virtual void jump_back( uintmax_t n ) { // jumps ahead the next n random numbers
		
		for ( size_t i = 0; i < N_SEEDS; i++ ) {
      
         Bitmatrix<uint32_t> A( MATRIX_INV[i] );
         _s[i] = ( A^n ) * _s[i];
      }
   }
   
   virtual void jump_back( uintmax_t e, uintmax_t c ) {   // jumps ahead the next n random numbers, where n = 2^e + c
   
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
   
   virtual void jump_cycle( void ) { // jump ahead an entire cycle of lfsr113
   
      const uint32_t A = 31, B = 29, C = 28, D = 25;
      jump_ahead( A + B + C + D, 0 );
      jump_back( A + B + C, 0 );
      jump_back( A + B + D, 0 );
      jump_back( A + C + D, 0 );
      jump_back( B + C + D, 0 );
      jump_ahead( A + B, 0 );
      jump_ahead( A + C, 0 );
      jump_ahead( A + D, 0 );
      jump_ahead( B + C, 0 );
      jump_ahead( B + D, 0 );
      jump_ahead( C + D, 0 );
      jump_back( A, 0 );
      jump_back( B, 0 );
      jump_back( C, 0 );
      jump_back( D, 0 );
      jump_ahead( 1 );
   }
/*
	virtual void jump_cycle( void ) { // jump ahead a full cycle of lfsr113
	
	   std::bitset<113> p( std::string( "11111111111111111111111101100110000000000000000000001101101011111111111111111111101001110000000000000000000000001" ) );
      for ( size_t i = 0; i < p.size(); ++i ) if ( p.test(i) ) jump_ahead( i, 0 );
	}
*/
   
   virtual uint32_t rng32( void ) {   // returns the next number (a 32-bit unsigned int)

      _s[0] = ( ( _s[0] & C1 ) << 18 ) ^ ( ( ( _s[0] <<  6 ) ^ _s[0] ) >> 13 );
      _s[1] = ( ( _s[1] & C2 ) <<  2 ) ^ ( ( ( _s[1] <<  2 ) ^ _s[1] ) >> 27 );
      _s[2] = ( ( _s[2] & C3 ) <<  7 ) ^ ( ( ( _s[2] << 13 ) ^ _s[2] ) >> 21 );
      _s[3] = ( ( _s[3] & C4 ) << 13 ) ^ ( ( ( _s[3] <<  3 ) ^ _s[3] ) >> 12 );
   
      return ( _s[0] ^ _s[1] ^ _s[2] ^ _s[3] ) & C0;
   }
   
   virtual uint64_t rng64( void ) {   // returns 64-bit integer
   
      uint64_t low  = rng32();
      uint64_t high = rng32();
      return low | ( high << 32 );
   }
   
   virtual double rng32_01( void ) {   // returns a double int the half-open interval [0,1)
   
      return double( rng32() ) * TWO32_INV;
   }
   
   virtual long double rng64_01( void ) {   // returns a long double in [0,1)
   
      return double( rng64() ) * TWO64_INV;
   }

private:

   uint32_t _s[ N_SEEDS ];

}; // end lfsr113 class
} // end namespace LFSR113

#endif // LFSR113_H
