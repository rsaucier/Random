// mod_math.h: modular math for 32-bit and 64-bit unsigned integers
// Ref: https://github.com/cmcqueen/simplerandom
// R. Saucier, 14 October 2016

#include <cstdint>
#include <cassert>

static const uint64_t    M         = 4294967296ULL;                     // 2^32
static const uint64_t    TWO32     = 4294967296ULL;                     // 2^32
static const long double TWO17     = 131072.0L;                         // 2^17
static const long double TWO35     = 34359738368.0L;                    // 2^35
static const long double TWO53     = 9007199254740992.0L;               // 2^53
static const long double TWO32_INV = 2.328306436538696289062500e-10L;   // 2^(-32)
static const long double TWO64_INV = 5.421010862427522170037264e-20L;   // 2^(-64)

// a + b mod m
uint32_t add_mod32( uint32_t a, uint32_t b, uint32_t m ) {

#ifdef UINT64_C // use the native 64-bit capability
   
   if ( b <= UINT32_MAX - a )
      return ( a + b ) % m;
   else
      return ( uint64_t( a ) + b ) % m;
   
#else // native 64-bit not available, so perform addition using 32-bit
   
   a %= m;
   b %= m;
   uint32_t t;
   if ( b <= UINT32_MAX - a )
      return ( a + b ) % m;
   
   if ( m <= ( UINT32_MAX >> 1 ) )
      return ( ( a % m ) + ( b % m ) ) % m;
   
   t = a + b;
   if ( t > uint32_t( m * 2 ) ) // m*2 must be truncated before compare
      t -= m;
   t -= m;
   return t % m;

#endif // UINT64_C
}

// a + b mod 2^32
uint32_t add32( uint32_t a, uint32_t b ) {

#ifdef UINT64_C // use the native 64-bit capability

   if ( b <= UINT32_MAX - a )
      return ( a + b ) % M;
   else
      return ( uint64_t( a ) + b ) % M;
   
#else // native 64-bit not available, so perform addition using 32-bit
   
   a %= M;
   b %= M;
   uint32_t t;
   if ( b <= UINT32_MAX - a )
      return ( a + b ) % M;
   
   if ( m <= ( UINT32_MAX >> 1 ) )
      return ( ( a % M ) + ( b % M ) ) % M;
   
   t = a + b;
   if ( t > uint32_t( M * 2 ) ) // m*2 must be truncated before compare
      t -= M;
   t -= M;
   return t % M;

#endif // UINT64_C
}

// a * b mod m
uint32_t mul_mod32( uint32_t a, uint32_t b, uint32_t m ) {

#ifdef UINT64_C // use the native 64-bit capability
   
   uint64_t t = uint64_t( a ) * b;
   return uint32_t( t % m );
   
#else // native 64-bit not available, so perform multiplication using 32-bit
   
   a %= m;
   b %= m;
   uint32_t r = 0;
   uint32_t t;
   
   if ( b >= m ) {
      
      if ( m > UINT32_MAX / 2u ) b -= m;
      else                       b %= m;
   }
   
   while ( a != 0 ) {
      
      if ( a & 1u ) {
         
         if ( b >= m - r ) r -= m;
         r += b;
      }
      a >>= 1u;
      
      t = b;
      if ( b >= m - t ) t -= m;
      b += t;
   }
   return r;
   
#endif // UINT64_C
}

// a * b mod 2^32
uint32_t mul32( uint32_t a, uint32_t b ) {

#ifdef UINT64_C // use the native 64-bit capability
   
   uint64_t t = uint64_t( a ) * b;
   return uint32_t( t % M );
   
#else // native 64-bit not available, so perform multiplication using 32-bit
   
   a %= M;
   b %= M;
   uint32_t r = 0;
   uint32_t t;
   
   if ( b >= m ) {
      
      if ( m > UINT32_MAX / 2u ) b -= M;
      else                       b %= M;
   }
   
   while ( a != 0 ) {
      
      if ( a & 1u ) {
         
         if ( b >= M - r ) r -= M;
         r += b;
      }
      a >>= 1u;
      
      t = b;
      if ( b >= M - t ) t -= M;
      b += t;
   }
   return r;
   
#endif // UINT64_C
}

// 32-bit methods

// a^n mod m
uint32_t pow_mod32( uint32_t a, uintmax_t n, uint32_t m ) {

   uint32_t r = 1;
   uint32_t t = a;
   
   for (;;) {
   
      if ( n & 1 ) r = mul_mod32( r, t, m );
      n >>= 1;
      if ( n == 0 ) break;
      t = mul_mod32( t, t, m );
   }
   return r;
}

// a^n mod m, where n = 2^e + c
uint32_t pow_mod32( uint32_t a, uintmax_t e, uintmax_t c, uint32_t m ) {
   
   if ( e == 0 ) return pow_mod32( a, c + 1, m );
   uint32_t t = a;
   for ( uintmax_t i = 0; i < e; ++i ) t = mul_mod32( t, t, m );
   return mul_mod32( pow_mod32( a, c, m ), t, m );
}

// a^n mod 2^32
uint32_t pow32( uint32_t a, uintmax_t n ) {

   uint32_t r = 1;
   uint32_t t = a;
   
   for (;;) {
   
      if ( n & 1 ) r *= t;
      n >>= 1;
      if ( n == 0 ) break;
      t *= t;
   }
   return r;
}

// a^n mod 2^32, where n = 2^e + c
uint64_t pow32( uint32_t a, uint32_t e, uint32_t c ) {
   
   if ( e == 0 ) return pow32( a, c + 1 );
   uint32_t t = a;
   for ( uint32_t i = 0; i < e; ++i ) t = mul32( t, t );
   return mul32( pow32( a, c ), t );
}

// sum first n terms of geometric series: 1 + a + ... + a^(n-1) mod m
uint32_t gs_mod32( uint32_t a, uint32_t n, uint32_t m ) {

   if ( n == 0 ) return 0;

   uint32_t t = a % m;
   uint32_t p = 1;
   uint32_t r = 0;

   while ( n > 1 ) {
   
      if ( n & 1 ) r = add_mod32( r, mul_mod32( p, pow_mod32( t, n - 1, m ), m ), m );
      p = mul_mod32( p, add_mod32( 1, t, m ), m );
      t = mul_mod32( t, t, m );
      n >>= 1;
   }
   r = add_mod32( r, p, m );
   return r;
}

// sum first n terms of geometric series: 1 + a + ... + a^(n-1) mod m, where n = 2^e + c
uint32_t gs_mod32( uint32_t a, uint32_t e, uint32_t c, uint32_t m ) {
   
   if ( e == 0 ) return gs_mod32( a, 1 + c, m );
   
   uint32_t t = a;
   uint32_t r = 1;

  for ( uint32_t i = 0; i < e; ++i ) {
      
      r = mul_mod32( r, add_mod32( 1, t, m ), m );
      t = mul_mod32( t, t, m );
   }
   if ( c == 0 ) return r;
   
   return add_mod32( r, mul_mod32( t, gs_mod32( a, c, m ), m ), m );
}

// sum first n terms of geometric series: 1 + a + ... + a^(n-1) mod 2^32
uint32_t gs32( uint32_t a, uintmax_t n ) {

   if ( n == 0 ) return 0;
   if ( n == 1 ) return 1;

   uint32_t t = a;
   uint32_t p = 1;
   uint32_t r = 0;

   while ( n > 1 ) {
   
      if ( n & 1 ) r += p * pow32( t, n - 1 );
      p *= ( 1 + t );
      t *= t;
      n >>= 1;
   }
   r += p;
   return r;
}

// sum first n terms of geometric series: 1 + a + ... + a^(n-1) mod 2^32, where n = 2^e + c
uint32_t gs32( uint32_t a, uint32_t e, uint32_t c ) {
   
   if ( e == 0 ) return gs32( a, 1 + c );
   
   uint32_t t = a;
   uint32_t r = 1;

  for ( uint32_t i = 0; i < e; ++i ) {
      
      r = mul32( r, add32( 1, t ) );
      t = mul32( t, t );
   }
   if ( c == 0 ) return r;
   
   return add32( r, mul32( t, gs32( a, c ) ) );
}

// 64-bit methods

#ifdef UINT64_C // the following require 64-bit capability

// 64-bit computation of a + b mod m
uint64_t add_mod64( uint64_t a, uint64_t b, uint64_t m ) {

   a %= m;
   b %= m;
   uint64_t t;
   if ( b <= UINT64_MAX - a )
      return ( a + b ) % m;

   if ( m <= ( UINT64_MAX >> 1 ) )
      return ( ( a % m ) + ( b % m ) ) % m;

   t = a + b;
   if ( t > uint64_t( m * 2 ) ) // m*2 must be truncated before compare
      t -= m;
   t -= m;
   return t % m;
}

// 64-bit computation of a + b mod 2^64
uint64_t add64( uint64_t a, uint64_t b ) {

   return a + b;
}

// 64-bit computation of a * b mod m
uint64_t mul_mod64( uint64_t a, uint64_t b, uint64_t m ) {

   uint64_t r = 0;
   uint64_t t;

   if ( b >= m ) {
   
      if ( m > UINT64_MAX / 2u ) b -= m;
      else                       b %= m;
   }

   while ( a != 0 ) {
   
      if ( a & 1 ) {
      
         if ( b >= m - r ) r -= m;
         r += b;
      }
      a >>= 1;

      t = b;
      if ( b >= m - t ) t -= m;
      b += t;
   }
   return r;
}

// 64-bit computation of a * b mod 2^64
uint64_t mul64( uint64_t a, uint64_t b ) {

   uint64_t r = 0;
   uint64_t t;

   while ( a != 0 ) {
   
      if ( a & 1 ) r += b;
      a >>= 1;
      t = b;
      b += t;
   }
   return r;
}

// 64-bit computation of a^n mod m
uint64_t pow_mod64( uint64_t a, uintmax_t n, uint64_t m ) {

   if ( n == 0 ) return 1;
   if ( n == 1 ) return a %= m;
   
   uint64_t r = 1;
   uint64_t t = a;
   
   for (;;) {
      
      if ( n & 1 ) r = mul_mod64( r, t, m );
      n >>= 1;
      if ( n == 0 ) break;
      t = mul_mod64( t, t, m );
   }
   return r;
}

// 64-bit computation of a^n mod m, where n = 2^e + c
uint64_t pow_mod64( uint64_t a, uintmax_t e, uintmax_t c, uint64_t m ) {
   
   if ( e == 0 ) return pow_mod64( a, c + 1, m );
   uint64_t t = a;
   for ( uint64_t i = 0; i < e; ++i ) t = mul_mod64( t, t, m );
   return mul_mod64( pow_mod64( a, c, m ), t, m );
}

// a^n mod 2^64
uint64_t pow64( uint64_t a, uintmax_t n ) {

   uint64_t r = 1;
   uint64_t t = a;
   
   for (;;) {
   
      if ( n & 1 ) r *= t;
      n >>= 1;
      if ( n == 0 ) break;
      t *= t;
   }
   return r;
}

// a^n mod 2^64, where n = 2^e + c
uint64_t pow64( uint64_t a, uint64_t e, uint64_t c ) {
   
   if ( e == 0 ) return pow64( a, c + 1 );
   uint64_t t = a;
   for ( uint64_t i = 0; i < e; ++i ) t = mul64( t, t );
   return mul64( pow64( a, c ), t );
}

// 64-bit sum first n terms of geometric series: 1 + a + ... + a^(n-1) mod m
uint64_t gs_mod64( uint64_t a, uintmax_t n, uint64_t m ) {

   if ( n == 0 ) return 0;

   uint64_t t = a % m;
   uint64_t p = 1;
   uint64_t r = 0;

   while ( n > 1 ) {
   
      if ( n & 1 ) r = add_mod64( r, mul_mod64( p, pow_mod64( t, n - 1, m ), m ), m );
      p = mul_mod64( p, add_mod64( 1, t, m ), m );
      t = mul_mod64( t, t, m );
      n >>= 1;
   }
   r = add_mod64( r, p, m );
   return r;
}

// 64-bit sum first n terms of geometric series: 1 + a + ... + a^(n-1) mod m, where n = 2^e + c
uint64_t gs_mod64( uint64_t a, uint32_t e, uint32_t c, uint64_t m ) {
   
   if ( e == 0 ) return gs_mod64( a, 1 + c, m );
   
   uint64_t t = a;
   uint64_t r = 1;

  for ( uint32_t i = 0; i < e; ++i ) {
      
      r = mul_mod64( r, add_mod64( 1, t, m ), m );
      t = mul_mod64( t, t, m );
   }
   if ( c == 0 ) return r;
   
   return add_mod64( r, mul_mod64( t, gs_mod64( a, c, m ), m ), m );
}

// 64-bit sum first n terms of geometric series: 1 + a + ... + a^(n-1) mod 2^64
uint64_t gs64( uint64_t a, uintmax_t n ) {

   if ( n == 0 ) return 0;
   if ( n == 1 ) return 1;

   uint64_t t = a;
   uint64_t p = 1;
   uint64_t r = 0;

   while ( n > 1 ) {
   
      if ( n & 1 ) r += mul64( p, pow64( t, n - 1 ) );
      p = mul64( p, 1 + t );
      t = mul64( t, t );
      n >>= 1;
   }
   r += p;
   return r;
}

// 64-bit sum first n terms of geometric series: 1 + a + ... + a^(n-1) mod 2^64, where n = 2^e + c
uint64_t gs64( uint64_t a, uint64_t e, uint64_t c ) {
   
   if ( e == 0 ) return gs64( a, 1 + c );
   
   uint64_t t = a;
   uint64_t r = 1;

  for ( uint32_t i = 0; i < e; ++i ) {
      
      r = mul64( r, add64( 1, t ) );
      t = mul64( t, t );
   }
   if ( c == 0 ) return r;
   
   return add64( r, mul64( t, gs64( a, c ) ) );
}

#endif // UINT64_C

// compute a + b mod m, where a, b and m must be < 2^35
double add_mod( double a, double b, double m ) {
   
   assert( a < TWO35 && b < TWO35 && m < TWO35 );
   double v = a + b;
   uintmax_t a1;
   
   if ( v >= TWO53 || v <= -TWO53 ) {
      a1 = static_cast<uintmax_t>( a / TWO17 );
      a -= a1 * TWO17;
      v  = a1;
      a1 = static_cast<uintmax_t>( v / m );
      v -= a1 * m;
      v = v * TWO17 + a + b;
   }
   
   a1 = static_cast<uintmax_t>( v / m );
   if ( ( v -= a1 * m ) < 0. ) return v += m;
   else                        return v;
}

// a * b mod m, where a, b, and m must be < 2^35
double mul_mod( double a, double b, double m ) {
   
   assert( a < TWO35 && b < TWO35 && m < TWO35 );
   double v = a * b;
   uintmax_t a1;
   
   if ( v >= TWO53 || v <= -TWO53 ) {
      a1 = static_cast<uintmax_t>( a / TWO17 );
      a -= a1 * TWO17;
      v  = a1 * b;
      a1 = static_cast<uintmax_t>( v / m );
      v -= a1 * m;
      v = v * TWO17 + a * b;
   }
   
   a1 = static_cast<uintmax_t>( v / m );
   if ( ( v -= a1 * m ) < 0. ) return v += m;
   else                        return v;
}

// compute a^n mod m
double pow_mod( double a, uintmax_t n, double m ) {
   
   if ( n == 0 ) return 1.;
   if ( n == 1 ) return a;
   
   double r = 1.;
   double t = a;

   for (;;) {
   
      if ( n & 1 ) r = mul_mod( r, t, m );
      n >>= 1;
      if ( n == 0 ) break;
      t = mul_mod( t, t, m );
   }
   return r;
}

// compute a^n mod m, where n = 2^e + c
uint64_t pow_mod( double a, uint32_t e, uintmax_t c, double m ) {
   
   if ( e == 0 ) return pow_mod( a, c + 1, m );
   double t = a;
   for ( uint32_t i = 0; i < e; ++i ) t = mul_mod( t, t, m );
   return mul_mod( pow_mod( a, c, m ), t, m );
}

// sum first n terms of geometric series: 1 + a + ... + a^(n-1) mod m
double gs_mod( double a, uintmax_t n, double m ) {
   
   if ( n == 0 ) return 0.;
   if ( n == 1 ) return 1.;
   
   double t = a;
   double p = 1.;
   double r = 0.;
   
   while ( n > 1 ) {
      
      if ( n & 1 ) r = add_mod( r, mul_mod( p, pow_mod( t, n - 1, m ), m ), m );
      p = mul_mod( p, 1. + t, m );
      t = mul_mod( t, t, m );
      n >>= 1;
   }
   r = add_mod( r, p, m );
   return r;
}

// sum first n terms of geometric series: 1 + a + ... + a^(n-1) mod m, where n = 2^e + c
double gs_mod( double a, uint32_t e, uint32_t c, double m ) {
      
   if ( e == 0 ) return gs_mod( a, 1 + c, m );
      
   double t = a;
   double r = 1;
      
   for ( uint32_t i = 0; i < e; ++i ) {
         
      r = mul_mod( r, add_mod( 1, t, m ), m );
      t = mul_mod( t, t, m );
   }
   if ( c == 0 ) return r;
      
   return add_mod( r, mul_mod( t, gs_mod( a, c, m ), m ), m );
}
