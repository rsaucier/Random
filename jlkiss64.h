// jlkiss64.h: Based upon Marsaglia's Keep It Simple Stupid RNG
// Ref: Good Practice in (Pseudo) Random Number Generation for Bioinformatics Applications
// David Jones, UCL Bioinformatics Group (d.jones@cs.ucl.ac.uk), May 7, 2010
// Period is (2^64)(2^64-1)(4294584393(2^31)-1)(698769069(2^31)-1) =
// 4,709,274,331,675,767,436,556,996,170,486,797,220,343,483,431,369,386,641,683,085,078,522,096,517,120
// or approximately 2^251 = 3.6 x 10^75
// Requires four 64-bit seeds to initialize.
// R. Saucier, 24 August 2016

#ifndef JLKISS64_H
#define JLKISS64_H

#include <iostream>
#include <iomanip>
using namespace std;

namespace JLKISS64 {
   
   static const bitmatrix64_t MATRIX = {
      {
         0x0008000440200011, 0x0010000880400022, 0x0020001100800044, 0x0040002201000088, 0x0080004402000110, 0x0100008804000220, 0x0200011008000440, 0x0400022010000880,
         0x0800044020001100, 0x1000088040002200, 0x2000110080004400, 0x4000220100008800, 0x8000440200011000, 0x0000880400022000, 0x0001100800044000, 0x0002201000088000,
         0x0004402000110000, 0x0008804040220001, 0x0011008080440002, 0x0022010100880004, 0x0044020201100008, 0x0088040402200010, 0x0110080804400020, 0x0220101008800040,
         0x0440202011000080, 0x0880404022000100, 0x1100808044000200, 0x2201010088000400, 0x4402020110000800, 0x8804040220001000, 0x1008080440002000, 0x2010100880004000,
         0x4020201100008000, 0x8040402200010000, 0x0080804400020000, 0x0101008800040000, 0x0202011000080000, 0x0404022000100000, 0x0808044000200000, 0x1010088000400000,
         0x2020110000800000, 0x4040220001000000, 0x8080440002000000, 0x0100080004000000, 0x0200100008000000, 0x0400200010000000, 0x0800400020000000, 0x1000800040000000,
         0x2001000080000000, 0x4002000100000000, 0x8004000200000000, 0x0008000400000000, 0x0010000800000000, 0x0020001000000000, 0x0040002000000000, 0x0080004000000000,
         0x0100008000000000, 0x0200010000000000, 0x0400020000000000, 0x0800040000000000, 0x1000080000000000, 0x2000100000000000, 0x4000200000000000, 0x8000400000000000
      }
   };
   static const bitmatrix64_t MATRIX_INV = {
      {
         0x90808c0404202201, 0x2101180808404402, 0x4202301010808804, 0x8404602021011008, 0x8880444402220011, 0x1100888804440022, 0x2201111008880044, 0x4402222011100088,
         0x8804444022200110, 0x1008888044400220, 0x2011110088800440, 0x4022220111000880, 0x8044440222001100, 0x0088880444002200, 0x0111100888004400, 0x0222201110008800,
         0x0444402220011000, 0x8888844440222001, 0x1111088880444002, 0x2222111100888004, 0x4444222201110008, 0x0888404402020011, 0x1110808804040022, 0x2221011008080044,
         0x4442022010100088, 0x8884044020200110, 0x1108088040400220, 0x2210110080800440, 0x4420220101000880, 0x8840440202001100, 0x1080880404002200, 0x2101100808004400,
         0x4202201010008800, 0x8404402020011000, 0x8880044400220001, 0x1100088800440002, 0x2200111000880004, 0x4400222001100008, 0x8800444002200010, 0x1000888004400020,
         0x2001110008800040, 0x4002220011000080, 0x8004440022000100, 0x0008880044000200, 0x0011100088000400, 0x0022200110000800, 0x0044400220001000, 0x0088800440002000,
         0x0111000880004000, 0x0222001100008000, 0x0444002200010000, 0x8888044400220001, 0x1110088800440002, 0x2220111000880004, 0x4440222001100008, 0x8880444002200010,
         0x1100888004400020, 0x2201110008800040, 0x4402220011000080, 0x8804440022000100, 0x1008880044000200, 0x2011100088000400, 0x4022200110000800, 0x8044400220001000
      }
   };
   static const uint64_t LC_MULT        = 0x14ada13ed78492ad;   // 1490024343005336237ULL;
   static const uint64_t LC_CONST       = 0x00000000075bcd15;   // 123456789ULL;
   static const uint64_t LC_MULT_INV    = 0xc5a2d1aa2af8a125;   // 14241175500494512421ULL;
   static const uint64_t MWC_MULT1      = 0x00000000fffa2849;   // 4294584393ULL;
   static const uint64_t MWC_MOD1       = 0xfffa2848ffffffff;   // 18445099517847011327ULL;
   static const uint64_t MWC_MULT1_INV  = 0x0000000100000000;   // 4294967296ULL;
   static const uint64_t MWC_MULT2      = 0x0000000029a65ead;   // 698769069ULL;
   static const uint64_t MWC_MOD2       = 0x29a65eacffffffff;   // 3001190298811367423ULL;
   static const uint64_t MWC_MULT2_INV  = 0x0000000100000000;   // 4294967296ULL;
   static const uint32_t N_SEEDS        = 4;                    // requires 4 64-bit words

class jlkiss64 : public Generator<uint64_t> {

public:
   jlkiss64( void ) { // default constructor
   }
   
   jlkiss64( std::vector<uint64_t> seed ) { // constructor from seed vector
   
      setState( seed );
   }
   
   virtual ~jlkiss64() { // default destructor
   
      std::cout << "deleting jlkiss64" << std::endl;
   }
   
   virtual void setState( std::vector<uint64_t> seed ) { // set the seeds from four 64-bit words
   
      assert( seed.size() >= N_SEEDS );
      _s1 = seed[0];
      _s2 = seed[1];
      _s3 = uint32_t( seed[2] >> 32 );
      _s4 = uint32_t( seed[2] );
      _s5 = uint32_t( seed[3] >> 32 );
      _s6 = uint32_t( seed[3] );
   }
   
   virtual void getState( std::vector<uint64_t>& seed ) { // get the seed vector
   
      assert( seed.size() >= N_SEEDS );
      seed[0] = _s1;
      seed[1] = _s2;
      seed[2] = ( uint64_t( _s3 ) << 32 ) + _s4;
      seed[3] = ( uint64_t( _s5 ) << 32 ) + _s6;
   }

   void jump_ahead( uintmax_t n ) { // jump ahead the next n random numbers
      
      _s1 = mul64( pow64( LC_MULT, n ), _s1 ) + mul64( LC_CONST, gs64( LC_MULT, n ) );
      
      Bitmatrix<uint64_t> A, B( MATRIX );
      A = B^n;
      _s2 = A * _s2;
      
      uint64_t a = _s3 + ( (uint64_t)_s4 << 32u );
      a = mul_mod64( pow_mod64( MWC_MULT1, n, MWC_MOD1 ), a, MWC_MOD1 );
      _s4 = ( uint32_t )( a >> 32u );
      _s3 = ( uint32_t )( a );
      
      a = _s5 + ( (uint64_t)_s6 << 32u );
      a = mul_mod64( pow_mod64( MWC_MULT2, n, MWC_MOD2 ), a, MWC_MOD2 );
      _s6 = ( uint32_t )( a >> 32u );
      _s5 = ( uint32_t )( a );
   }

   void jump_ahead( uintmax_t e, uintmax_t c ) {   // jump ahead the next n random numbers, where n = 2^e + c
      
      _s1 = mul64( pow64( LC_MULT, e, c ), _s1 ) + mul64( LC_CONST, gs64( LC_MULT, e, c ) );
      
      Bitmatrix<uint64_t> A, B;
   
      if ( e ) {
         B = MATRIX;
         for ( uint64_t i = 0; i < e; i++ ) B *= B;
      }
      A = MATRIX;
      A = A^c;
      if ( e ) A *= B;
      _s2 = A * _s2;
      
      uint64_t a = _s3 + ( (uint64_t)_s4 << 32u );
      a = mul_mod64( pow_mod64( MWC_MULT1, e, c, MWC_MOD1 ), a, MWC_MOD1 );
      _s4 = ( uint32_t )( a >> 32u );
      _s3 = ( uint32_t )( a );
      
      a = _s5 + ( (uint64_t)_s6 << 32u );
      a = mul_mod64( pow_mod64( MWC_MULT2, e, c, MWC_MOD2 ), a, MWC_MOD2 );
      _s6 = ( uint32_t )( a >> 32u );
      _s5 = ( uint32_t )( a );
   }
   
   void jump_back( uintmax_t n ) { // jump back the next n random numbers
   
      _s1 = mul64( pow64( LC_MULT_INV, n ), _s1 - LC_CONST ) + LC_CONST - mul64( LC_CONST, gs64( LC_MULT_INV, n ) );
      
      Bitmatrix<uint64_t> A, B( MATRIX_INV );
      A = B^n;
      _s2 = A * _s2;
      
      uint64_t a = _s3 + ( (uint64_t)_s4 << 32u );
      a = mul_mod64( pow_mod64( MWC_MULT1_INV, n, MWC_MOD1 ), a, MWC_MOD1 );
      _s4 = ( uint32_t )( a >> 32u );
      _s3 = ( uint32_t )( a );
      
      a = _s5 + ( (uint64_t)_s6 << 32u );
      a = mul_mod64( pow_mod64( MWC_MULT2_INV, n, MWC_MOD2 ), a, MWC_MOD2 );
      _s6 = ( uint32_t )( a >> 32u );
      _s5 = ( uint32_t )( a );
   }
   
   void jump_back( uintmax_t e, uintmax_t c ) {   // jump ahead the next n random numbers, where n = 2^e + c
   
      _s1 = mul64( pow64( LC_MULT_INV, e, c ), _s1 - LC_CONST ) + LC_CONST - mul64( LC_CONST, gs64( LC_MULT_INV, e, c ) );
      
      Bitmatrix<uint64_t> A, B;
   
      if ( e ) {
         B = MATRIX_INV;
         for ( uint64_t i = 0; i < e; i++ ) B *= B;
      }
      A = MATRIX_INV;
      A = A^c;
      if ( e ) A *= B;
      _s2 = A * _s2;
      
      uint64_t a = _s3 + ( (uint64_t)_s4 << 32u );
      a = mul_mod64( pow_mod64( MWC_MULT1_INV, e, c, MWC_MOD1 ), a, MWC_MOD1 );
      _s4 = ( uint32_t )( a >> 32u );
      _s3 = ( uint32_t )( a );
      
      a = _s5 + ( (uint64_t)_s6 << 32u );
      a = mul_mod64( pow_mod64( MWC_MULT2_INV, e, c, MWC_MOD2 ), a, MWC_MOD2 );
      _s6 = ( uint32_t )( a >> 32u );
      _s5 = ( uint32_t )( a );
   }
   
   virtual void jump_cycle( void ) { // jump ahead an entire cycle of JLKISS64

      jump_ahead( 251, 0 );
      jump_ahead( 249, 0 );
      jump_ahead( 246, 0 );
      jump_ahead( 245, 0 );
      jump_ahead( 243, 0 );
      jump_ahead( 240, 0 );
      jump_ahead( 238, 0 );
      jump_ahead( 236, 0 );
      jump_ahead( 235, 0 );
      jump_ahead( 233, 0 );
      jump_ahead( 231, 0 );
      jump_ahead( 230, 0 );
      jump_ahead( 228, 0 );
      jump_ahead( 226, 0 );
      jump_ahead( 224, 0 );
      jump_ahead( 221, 0 );
      jump_ahead( 219, 0 );
      jump_ahead( 216, 0 );
      jump_ahead( 215, 0 );
      jump_ahead( 214, 0 );
      jump_ahead( 211, 0 );
      jump_ahead( 209, 0 );
      jump_ahead( 208, 0 );
      jump_ahead( 207, 0 );
      jump_ahead( 200, 0 );
      jump_ahead( 199, 0 );
      jump_ahead( 198, 0 );
      jump_ahead( 196, 0 );
      jump_ahead( 194, 0 );
      jump_ahead( 191, 0 );
      jump_ahead( 189, 0 );
      jump_ahead( 183, 0 );
      jump_ahead( 182, 0 );
      jump_ahead( 178, 0 );
      jump_ahead( 177, 0 );
      jump_ahead( 174, 0 );
      jump_ahead( 173, 0 );
      jump_ahead( 168, 0 );
      jump_ahead( 167, 0 );
      jump_ahead( 165, 0 );
      jump_ahead( 163, 0 );
      jump_ahead( 162, 0 );
      jump_ahead( 161, 0 );
      jump_ahead( 160, 0 );
      jump_ahead( 159, 0 );
      jump_ahead( 158, 0 );
      jump_ahead( 156, 0 );
      jump_ahead( 154, 0 );
      jump_ahead( 153, 0 );
      jump_ahead( 149, 0 );
      jump_ahead( 148, 0 );
      jump_ahead( 146, 0 );
      jump_ahead( 142, 0 );
      jump_ahead( 141, 0 );
      jump_ahead( 140, 0 );
      jump_ahead( 139, 0 );
      jump_ahead( 138, 0 );
      jump_ahead( 137, 0 );
      jump_ahead( 133, 0 );
      jump_ahead( 131, 0 );
      jump_ahead( 130, 0 );
      jump_ahead( 126, 0 );
      jump_ahead( 124, 0 );
      jump_ahead( 122, 0 );
      jump_ahead( 119, 0 );
      jump_ahead( 118, 0 );
      jump_ahead( 116, 0 );
      jump_ahead( 110, 0 );
      jump_ahead( 105, 0 );
      jump_ahead( 104, 0 );
      jump_ahead( 102, 0 );
      jump_ahead( 101, 0 );
      jump_ahead( 100, 0 );
      jump_ahead( 99, 0 );
      jump_ahead( 97, 0 );
      jump_ahead( 95, 0 );
      jump_ahead( 94, 0 );
      jump_ahead( 93, 0 );
      jump_ahead( 92, 0 );
      jump_ahead( 91, 0 );
      jump_ahead( 90, 0 );
      jump_ahead( 89, 0 );
      jump_ahead( 88, 0 );
      jump_ahead( 87, 0 );
      jump_ahead( 86, 0 );
      jump_ahead( 85, 0 );
      jump_ahead( 84, 0 );
      jump_ahead( 83, 0 );
      jump_ahead( 82, 0 );
      jump_ahead( 81, 0 );
      jump_ahead( 80, 0 );
      jump_ahead( 79, 0 );
      jump_ahead( 78, 0 );
      jump_ahead( 77, 0 );
      jump_ahead( 76, 0 );
      jump_ahead( 75, 0 );
      jump_ahead( 74, 0 );
      jump_ahead( 73, 0 );
      jump_ahead( 72, 0 );
      jump_ahead( 71, 0 );
      jump_ahead( 70, 0 );
      jump_ahead( 69, 0 );
      jump_ahead( 68, 0 );
      jump_ahead( 67, 0 );
      jump_ahead( 66, 0 );
      jump_ahead( 65, 0 );
      jump_ahead( 64, 0 );
   }
/*
	virtual void jump_cycle( void ) { // jump ahead a full cycle of jlkiss
	
	   std::bitset<251> p( std::string( "101001101001010110101101010100101001110010111000000111010100101000001100011001100001101011111101011000110100011111100010110001010100110100000100001101111010111111111111111111111111111111110000000000000000000000000000000000000000000000000000000000000000" ) );
      for ( size_t i = 0; i < p.size(); ++i ) if ( p.test(i) ) jump_ahead( i, 0 );
	}
*/

   uint32_t rng32( void ) { // returns the next random number (as a 32-bit unsigned int)

      _s1 = LC_MULT * _s1 + LC_CONST;

      _s2 ^= ( _s2 << 21 ), _s2 ^= ( _s2 >> 17 ), _s2 ^= ( _s2 << 30 );

      uint64_t a = MWC_MULT1 * _s3 + _s4;
      _s4 = ( a >> 32u );
      _s3 = uint32_t( a );

      //return (uint32_t)( _s1 >> 32 ) + (uint32_t)_s2 + _s3;
      return ( uint32_t )( _s1 + _s2 + _s3 );
   }

   uint64_t rng64( void ) { // returns the next random number (as a 64-bit unsigned integer)

      _s1 = LC_MULT * _s1 + LC_CONST;

      _s2 ^= ( _s2 << 21 ), _s2 ^= ( _s2 >> 17 ), _s2 ^= ( _s2 << 30 );

      uint64_t a = MWC_MULT1 * _s3 + _s4;
      _s4 = ( a >> 32u );    // upper 32 bits of a
      _s3 = uint32_t( a );   // lower 32 bits of a
      
      a = MWC_MULT2 * _s5 + _s6;
      _s6 = ( a >> 32u );    // upper 32 bits of a
      _s5 = uint32_t( a );   // lower 32 bits of a

      // use _s3 for lower 32 bits and _s5 << 32 for upper 32 bits of 64-bit word
      return _s1 + _s2 + _s3 + ( ( uint64_t )_s5 << 32 );
   }

   double rng32_01( void ) { // returns a random number in the half-open interval [0,1)
   
      return double( rng32() ) * TWO32_INV;
   }

   long double rng64_01( void ) { // returns a random number in the half-open interval [0,1)
   
      return ( long double )( rng64() ) * TWO64_INV;
   }
/*
uint32_t myrng32( void ) {

   uint64_t _s1 = 123456789123ULL;
   uint64_t _s2 = 987654321987ULL;
   uint32_t _s3 = 43219876, _s4 = 6543217;
   
   uint64_t t;
   
   _s1 = mul64( 1490024343005336237ULL, _s1 ) + 123456789;
   
   _s2 ^= ( _s2 << 21 ), _s2 ^= ( _s2 >> 17 ), _s2 ^= ( _s2 << 30 );
   
   t = mul64( 4294584393ULL, _s3 ) + _s4;
   _s4 = t >> 32;
   _s3 = t;
   
   return ( uint32_t )( _s1 >> 32 )  + ( uint32_t )_s2 + _s3;
}

uint64_t myrng64( void ) {

   uint64_t _s1 = 123456789123ULL;
   uint64_t _s2 = 987654321987ULL;
   uint32_t _s3 = 43219876, _s4 = 6543217, _s5 = 21987643, _s6 = 1732654;
   
   uint64_t t;
   
   _s1 = mul64( 1490024343005336237ULL, _s1 ) + 123456789;
   
   _s2 ^= ( _s2 << 21 ), _s2 ^= ( _s2 >> 17 ), _s2 ^= ( _s2 << 30 );
   
   t = 4294584393ULL * _s3 + _s4;
   _s4 = t >> 32;
   _s3 = t;
   
   t = 4246477509ULL * s5 + s6;
   s6 = t >> 32;
   s5 = t;
   
   return _s1 + _s2 + _s3 + ( ( uint64_t )s5 << 32 );
}
*/

private:
   
   uint64_t _s1, _s2;             // two 64-bit
   uint32_t _s3, _s4, _s5, _s6;   // important that these be 32-bit and not 64-bit

}; // end jlkiss64 class
} // end namespace JLKISS64

#endif // JLKISS64_H
