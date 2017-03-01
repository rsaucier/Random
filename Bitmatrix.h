// Bitmatrix.h: template class for 32 x 32 or 64 x 64 matrices using mod 2 arithmetic
// R. Saucier, August 2016

#ifndef BITMATRIX_H
#define BITMATRIX_H

#include <cstdint>   // for uint32_t and uint64_t
#include <climits>   // for CHAR_BIT, the number of bits per byte

typedef struct { uint32_t matrix[32]; } bitmatrix32_t;
typedef struct { uint64_t matrix[64]; } bitmatrix64_t;

template <class T>
class Bitmatrix {

public:

   static const unsigned int N_BITS = CHAR_BIT * sizeof( T );   // number of bits

public:

   Bitmatrix( void ) {   // default constructor
   }
   
   Bitmatrix( const bitmatrix32_t& A ) {   // constructor from array of 32-bit constants
   
      for ( T i = 0; i < N_BITS; i++ ) _matrix[i] = A.matrix[i];
   }
   
   Bitmatrix( const bitmatrix64_t& A ) {   // constructor from array of 64-bit constants
   
      for ( T i = 0; i < N_BITS; i++ ) _matrix[i] = A.matrix[i];
   }
   
  ~Bitmatrix( void ) {   // default destructor
   }
  
   Bitmatrix( const Bitmatrix& A ) {   // copy constructor
   
      for ( T i = 0; i < N_BITS; i++ ) _matrix[i] = A._matrix[i];
   }
   
   Bitmatrix& operator=( const Bitmatrix& A ) {   // assignment operator
   
      if ( this != &A ) for ( T i = 0; i < N_BITS; i++ ) _matrix[i] = A._matrix[i];
      return *this;
   }
   
   void identity( Bitmatrix& A ) {   // create an identity matrix
      
      T v = T(1);
      for ( T i = 0; i < N_BITS; i++, v <<= 1 ) A._matrix[i] = v;
   }
    
   T matrix( T i ) {   // return the ith vector of the bitmatrix
   
      return _matrix[i];
   }
   
   // overloaded operators

   friend T operator*( const Bitmatrix<T>& A, T v ) {   // matrix multiplication of a vector
   
      T r = T(0);
      T b = T(1);
      for ( T i = 0; i < N_BITS; i++, v >>= 1 ) if ( v & b ) r ^= A._matrix[i];
      return r;
   }
   
   friend Bitmatrix operator*( Bitmatrix<T>& A, const Bitmatrix<T>& B ) {   // multiplication of two Bitmatrices
   
      Bitmatrix<T> C;
   
      for ( T i = 0; i < N_BITS; i++ ) C._matrix[i] = A * B._matrix[i];
      return C;
   }
   
   Bitmatrix& operator*=( const Bitmatrix<T>& A ) {   // multiplication assignment of two Bitmatrices
   
      return *this = *this * A;
   }
   
   Bitmatrix operator^( T n ) {   // return A^n, Bitmatrix A to the power n
   
      Bitmatrix<T> B, A = *this;
   
      identity( B );
      T b = T(1);
   
      while ( n > 0 ) {   // Knuth's "exponentiation by squaring" algorithm
      
         if ( n & b ) B *= A;
         A *= A;
         n >>= 1;
      }
      return B;
   }
   
   friend Bitmatrix pow( Bitmatrix<T>& A, T e, T c ) {   // return A^n, Bitmatrix A to the power n, where n = 2^e + c
      
      Bitmatrix<T> B;//, A = *this;
      if ( e > 0 ) {
         B = A;
         for ( T i = 0; i < e; i++ ) B *= B;
      }
      A = A^c;
      if ( e ) A *= B;
      return A;
   }
/*
   // multiply matrix A times vector v and return the vector result, b = A * v
   friend T bitmatrix_mul( const Bitmatrix<T>& A, T v ) {
      
      T b = T(0);
      for ( T i = 0; i < N_BITS; i++, v >>= 1 ) if ( v & 1 ) b ^= A._matrix[i];
      return b;
   }

   // multiply two matrices, A * B and store in A, so that A *= B
   friend void bitmatrix_mul( Bitmatrix<T>& A, const Bitmatrix<T>& B ) {
   
      Bitmatrix<T> C;
   
      for ( T i = 0; i < N_BITS; i++ ) C._matrix[i] = A * B._matrix[i];
      for ( T i = 0; i < N_BITS; i++ ) A._matrix[i] = C._matrix[i];
   }

   // raise matrix A to the power n and store in B matrix, so that B = A^n
   // Knuth's algorithm is "exponentiation by squaring" and time complexity is O(log n)
   friend void bitmatrix_pow( Bitmatrix<T>& B, const Bitmatrix<T>& A, uintmax_t n ) {
   
      Bitmatrix<T> C;
   
      for ( T i = 0; i < N_BITS; i++ ) C._matrix[i] = A._matrix[i];
      identity( B );
   
      while ( n > 0 ) {
      
         if ( n & 1 ) B *= C;
         C *= C;
         n >>= 1;
      }
   }
*/
private:

   T _matrix[N_BITS];
};
   // declaration of friends
   //void identity( Bitmatrix<uint32_t>& A );
   uint32_t operator*( const Bitmatrix<uint32_t>& A, uint32_t v );
   Bitmatrix<uint32_t> operator*( Bitmatrix<uint32_t>& A, const Bitmatrix<uint32_t>& B );
   Bitmatrix<uint32_t> pow( Bitmatrix<uint32_t>& A, uint32_t e, uint32_t c );

   //void identity( Bitmatrix<uint64_t>& A );
   uint64_t operator*( const Bitmatrix<uint64_t>& A, uint32_t v );
   Bitmatrix<uint64_t> operator*( Bitmatrix<uint64_t>& A, const Bitmatrix<uint64_t>& B );


#endif // BITMATRIX_H
