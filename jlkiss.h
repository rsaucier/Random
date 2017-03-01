// jlkiss.h: Based upon Marsaglia's Keep It Simple Stupid RNG
// Ref: Good Practice in (Pseudo) Random Number Generation for Bioinformatics Applications
// David Jones, UCL Bioinformatics Group (d.jones@cs.ucl.ac.uk), May 7, 2010
// Cycle length is (2^64)(2^64-1)(4294584393(2^31)-1) = 3138271061012620924047441856806230331094853687768430673920,
// or approximately 2^191
// Period of MWC is 4294584393(2^31)-1 = 9222549758923505663, where 4294584393 = 2^32 - 2^19 + 2^17 + 2^13 + 2^11 + 2^6 + 2^3 + 1
// R. Saucier, 30 August 2016

#ifndef JLKISS_H
#define JLKISS_H

namespace JLKISS {
   
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
   static const uint64_t LC_MULT      = 0x14ada13ed78492adull;   // 1490024343005336237ull;
   static const uint64_t LC_CONST     = 0x00000000075bcd15ull;   // 123456789ull;
   static const uint64_t LC_MULT_INV  = 0xc5a2d1aa2af8a125ull;   // 14241175500494512421ull;
   static const uint64_t MWC_MULT     = 0x00000000fffa2849ull;   // 4294584393ull;
   static const uint64_t MWC_MOD      = 0xfffa2848ffffffffull;   // 18445099517847011327ull;
   static const uint64_t MWC_MULT_INV = 0x0000000100000000ull;   // 4294967296ull;
   static const uint64_t SR_PERIOD    = 0xffffffffffffffffull;   // 18446744073709551615ull;
   static const uint64_t MWC_PERIOD   = 0x7ffd14247fffffffull;   // 9222549758923505663ull;
   static const uint32_t N_SEEDS      = 3;                       // requires three 64-bit words

class jlkiss : public Generator<uint64_t> {

public:
   jlkiss( void ) { // default constructor
   }
   
   jlkiss( std::vector<uint64_t> seed ) { // constructor from seed vector
   
      setState( seed );
   }
   
   virtual ~jlkiss() { // default destructor
   
      std::cout << "deleting jlkiss" << std::endl;
   }
   
   virtual void setState( std::vector<uint64_t> seed ) { // set the seeds
   
      assert( seed.size() >= N_SEEDS );
      _s1 = seed[0];
      _s2 = seed[1];
      _s3 = uint32_t( seed[2] >> 32 );
      _s4 = uint32_t( seed[2] );
   }
   
   virtual void getState( std::vector<uint64_t>& seed ) { // get the seed vector
   
      assert( seed.size() >= N_SEEDS );
      seed[0] = _s1;
      seed[1] = _s2;
      seed[2] = ( uint64_t( _s3 ) << 32 ) + _s4;
   }

   void jump_ahead( uintmax_t n ) { // jump ahead the next n random numbers
      
      _s1 = mul64( pow64( LC_MULT, n ), _s1 ) + mul64( LC_CONST, gs64( LC_MULT, n ) );
      
      Bitmatrix<uint64_t> A, B( MATRIX );
      A = B^n;
      _s2 = A * _s2;
      
      uint64_t a = _s3 + ( (uint64_t)_s4 << 32u );
      a = mul_mod64( pow_mod64( MWC_MULT, n, MWC_MOD ), a, MWC_MOD );
      _s4 = ( uint32_t )( a >> 32u );
      _s3 = ( uint32_t )( a );
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
      a = mul_mod64( pow_mod64( MWC_MULT, e, c, MWC_MOD ), a, MWC_MOD );
      _s4 = ( uint32_t )( a >> 32u );
      _s3 = ( uint32_t )( a );
   }
   
   void jump_back( uintmax_t n ) { // jump back the next n random numbers
   
      _s1 = mul64( pow64( LC_MULT_INV, n ), _s1 - LC_CONST ) + LC_CONST - mul64( LC_CONST, gs64( LC_MULT_INV, n ) );
      
      Bitmatrix<uint64_t> A, B( MATRIX_INV );
      A = B^n;
      _s2 = A * _s2;
      
      uint64_t a = _s3 + ( (uint64_t)_s4 << 32u );
      a = mul_mod64( pow_mod64( MWC_MULT_INV, n, MWC_MOD ), a, MWC_MOD );
      _s4 = ( uint32_t )( a >> 32u );
      _s3 = ( uint32_t )( a );
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
      a = mul_mod64( pow_mod64( MWC_MULT_INV, e, c, MWC_MOD ), a, MWC_MOD );
      _s4 = ( uint32_t )( a >> 32u );
      _s3 = ( uint32_t )( a );
   }
   
   virtual void jump_cycle( void ) { // jump ahead an entire cycle of JLKISS

      jump_ahead( 191, 0 ); jump_back( 127, 0 );
      jump_back( 178, 0 );  jump_ahead( 114, 0 );
      jump_ahead( 176, 0 ); jump_back( 112, 0 );
      jump_ahead( 172, 0 ); jump_back( 108, 0 );
      jump_ahead( 170, 0 ); jump_back( 106, 0 );
      jump_ahead( 165, 0 ); jump_back( 101, 0 );
      jump_ahead( 162, 0 ); jump_back( 98, 0 );
      jump_ahead( 159, 0 ); jump_back( 95, 0 );
      jump_back( 128, 0 );  jump_ahead( 64, 0 );
   }
/*
	virtual void jump_cycle( void ) { // jump ahead a full cycle of jlkiss
	
	   std::bitset<191> p( std::string( "11111111111110100010100001001000111111111111111111111111111111010000000000000101110101111011011100000000000000000000000000000010000000000000000000000000000000000000000000000000000000000000000" ) );
      for ( size_t i = 0; i < p.size(); ++i ) if ( p.test(i) ) jump_ahead( i, 0 );
	}
*/

   uint32_t rng32( void ) {

      _s1 = LC_MULT * _s1 + LC_CONST;

      _s2 ^= ( _s2 << 21 ), _s2 ^= ( _s2 >> 17 ), _s2 ^= ( _s2 << 30 );

      uint64_t a = MWC_MULT * _s3 + _s4;
      _s4 = ( a >> 32u );
      _s3 = uint32_t( a );

      //return (uint32_t)_s1 + (uint32_t)_s2 + _s3;   // this didn't work
      //return (uint32_t)( _s1 >> 32 ) + (uint32_t)_s2 + _s3;   // this doesn't seem to work
		return _s1 + _s2 + _s3;   // this passes all the tests, small crush, crush, and big crush
   }

   uint64_t rng64( void ) {

      _s1 = LC_MULT * _s1 + LC_CONST;

      _s2 ^= ( _s2 << 21 ), _s2 ^= ( _s2 >> 17 ), _s2 ^= ( _s2 << 30 );

      uint64_t a = MWC_MULT * _s3 + _s4;
      _s4 = ( a >> 32u );
      _s3 = uint32_t( a );

      return _s1 + _s2 + ( (uint64_t)_s4 << 32 ) + _s3;
   }

   double rng32_01( void ) { // returns a random number in the half-open interval [0,1)
   
      return double( rng32() ) * TWO32_INV;
   }

   long double rng64_01( void ) { // returns a random number in the half-open interval [0,1)
   
      return ( long double )( rng64() ) * TWO64_INV;
   }

private:
   
   uint64_t _s1, _s2;
   uint32_t _s3, _s4;

}; // end jlkiss class
} // end namespace JLKISS

#endif // JLKISS_H
