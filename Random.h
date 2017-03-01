// Random.h: Definition and Implementation of Random Number Distribution Class
// This rewrite of the following reference decouples the distributions from the generators
// Ref: Richard Saucier, "Computer Generation of Statistical Distributions," ARL-TR-2168,
//      US Army Research Laboratory, Aberdeen Proving Ground, MD, 21005-5068, March 2000.

#ifndef RANDOM_H
#define RANDOM_H

#include "Generator.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <functional>
#include <cassert>
#include <cmath>
#include <climits>
#include <cfloat>     // for FLT_MIN and FLT_MAX
#include <unistd.h>   // for getpid
#include <map>

namespace rnd {   // rnd namespace

// for convenience, define some data structures for bivariate distributions

struct point2d {   // cartesian coordinates in 2-D

   double x, y;
   point2d& operator+=( const point2d& p ) {
      x += p.x;
      y += p.y;
      return *this;
   }
   point2d& operator-=( const point2d& p ) {
      x -= p.x;
      y -= p.y;
      return *this;
   }
   point2d& operator*=( double scalar ) {
      x *= scalar;
      y *= scalar;
      return *this;
   }
   point2d& operator/=( double scalar ) {
      x /= scalar;
      y /= scalar;
      return *this;
   }
};

struct point3d {  // spherical coordinates on unit sphere

   point3d( void ) {
   }
   
   point3d( double th, double ph ) : theta( th ), phi( ph ) {
   }

  ~point3d( void ) {
   }
   
   double theta, phi;
   double x( void ) { return sin( theta ) * cos( phi ); }   // x-coordinate
   double y( void ) { return sin( theta ) * sin( phi ); }   // y-coordinate
   double z( void ) { return cos( theta ); }                // z-coordinate
};

template <class Typename>   // for 32-bit and 64-bit generators
class Random {

public:

   Random( Generator<Typename> *gen ) { _gen = gen; }
  ~Random( void ) {}  // default destructor

// Continuous Distributions

   double arcsine( double xMin = 0., double xMax = 1. ) { // Arc Sine

      double q = sin( M_PI_2 * _u() );
      return xMin + ( xMax - xMin ) * q * q;
   }
   
   double beta( double v, double w,                    // Beta
                double xMin = 0., double xMax = 1. ) { // (v > 0. and w > 0.)

      if ( v < w ) return xMax - ( xMax - xMin ) * beta( w, v );
      double y1 = gamma( 0., 1., v );
      double y2 = gamma( 0., 1., w );
      return xMin + ( xMax - xMin ) * y1 / ( y1 + y2 );
   }
   
   double cauchy( double a = 0., double b = 1. ) { // Cauchy (or Lorentz)

      // a is the location parameter and b is the scale parameter
      // b is the half width at half maximum (HWHM) and variance doesn't exist
   
      assert( b > 0. );
      
      return a + b * tan( M_PI * uniform( -0.5, 0.5 ) );
   }

   double chiSquare( int df ) { // Chi-Square
      
      assert( df >= 1 );
      
      return gamma( 0., 2., 0.5 * double( df ) );
   }

   double cosine( double xMin = 0., double xMax = 1. ) { // Cosine
   
      assert( xMin < xMax );
   
      double a = 0.5 * ( xMin + xMax );    // location parameter
      double b = ( xMax - xMin ) / M_PI;   // scale parameter
   
      return a + b * asin( uniform( -1., 1. ) );
   }
   
   double doubleLog( double xMin = -1., double xMax = 1. ) { // Double Log

      assert( xMin < xMax );
   
      double a = 0.5 * ( xMin + xMax );    // location parameter
      double b = 0.5 * ( xMax - xMin );    // scale parameter
   
      if ( bernoulli( 0.5 ) ) return a + b * _u() * _u();
      else                    return a - b * _u() * _u();
   }

   double erlang( double b, int c ) { // Erlang (b > 0. and c >= 1)

      assert( b > 0. && c >= 1 );
  
      double prod = 1.;
      for ( int i = 0; i < c; i++ ) prod *= _u();
      
      return -b * log( prod );
   }
   
   double exponential( double a = 0., double c = 1. ) { // Exponential
                                                    // location a, shape c
      assert( c > 0.0 );
   
      return a - c * log( _u() );
   }
   
   double extremeValue( double a = 0., double c = 1. ) { // Extreme Value
                                                     // location a, shape c
      assert( c > 0. );
   
      return a + c * log( -log( _u() ) );
   }
   
   double fRatio( int v, int w ) { // F Ratio (v and w >= 1)

      assert( v >= 1 && w >= 1 );
   
      return ( chiSquare( v ) / v ) / ( chiSquare( w ) / w );
   }
   
   double gamma( double a, double b, double c ) { // Gamma
                                                 // location a, scale b, shape c
      assert( b > 0. && c > 0. );
   
      if ( c < 1. ) {
         const double C = 1. + c / M_E;
         while ( true ) {
            double p = C * _u();      
            if ( p > 1. ) {
               double y = -log( ( C - p ) / c );
               if ( _u() <= pow( y, c - 1. ) ) return a + b * y;
            }
            else {
               double y = pow( p, 1. / c );
               if ( _u() <= exp( -y ) ) return a + b * y;
            }
         }
      }
      else if ( c == 1.0 ) return exponential( a, b );
      else {
         const double A = 1. / sqrt( 2. * c - 1. );
         const double B = c - log( 4. );
         const double Q = c + 1. / A;
         const double T = 4.5;
         const double D = 1. + log( T );
         while ( true ) {
            double p1 = _u();
            double p2 = _u();
            double v = A * log( p1 / ( 1. - p1 ) );
            double y = c * exp( v );
            double z = p1 * p1 * p2;
            double w = B + Q * v - y;
            if ( w + D - T * z >= 0. || w >= log( z ) ) return a + b * y;
         }
      }
   }
   
   double laplace( double a = 0., double b = 1. ) { // Laplace
                                                // (or double exponential)
      assert( b > 0. );

      // composition method
  
      if ( bernoulli( 0.5 ) ) return a + b * log( _u() );
      else                    return a - b * log( _u() );
   }
   
   double logarithmic( double xMin = 0., double xMax = 1. ) { // Logarithmic

      assert( xMin < xMax );
   
      double a = xMin;          // location parameter
      double b = xMax - xMin;   // scale parameter
   
      // use convolution formula for product of two IID uniform variates

      return a + b * _u() * _u();
   }
   
   double logistic( double a = 0., double c = 1. ) { // Logistic

      assert( c > 0. );

      return a - c * log( 1. / _u() - 1. );
   }
   
   double lognormal( double a, double mu, double sigma ) { // Lognormal

      return a + exp( normal( mu, sigma ) );
   }
   
   double normal( double mu = 0., double sigma = 1. ) { // Normal

      assert( sigma > 0. );
   
      static bool f = true;
      static double p2, q;
      double p1, p;
   
      if ( f ) {
         do {
            p1 = uniform( -1., 1. );
            p2 = uniform( -1., 1. );
            p = p1 * p1 + p2 * p2;
         } while ( p >= 1. );
         f = false;
         q = sqrt( -2. * log( p ) / p );
         return mu + sigma * p1 * q;
      }
      f = true;
      return mu + sigma * p2 * q;
   }
   
   double parabolic( double xMin = 0., double xMax = 1. ) { // Parabolic
  
      assert( xMin < xMax );
   
      double a    = 0.5 * ( xMin + xMax );        // location parameter
      double yMax = _parabola( a, xMin, xMax );   // maximum function range
   
      return userSpecified( _parabola, xMin, xMax, 0., yMax );
   }
   
   double pareto( double c ) { // Pareto
                               // shape c
      assert( c > 0. );
   
      return pow( _u(), -1. / c );
   }
   
   double pearson5( double b, double c ) { // Pearson Type 5
                                           // scale b, shape c
      assert( b > 0. && c > 0. );
   
      return 1. / gamma( 0., 1. / b, c );
   }
   
   double pearson6( double b, double v, double w ) { // Pearson Type 6
                                                     // scale b, shape v & w
      assert( b > 0. && v > 0. && w > 0. );
   
      return gamma( 0., b, v ) / gamma( 0., b, w );
   }
   
   double power( double c ) { // Power
                              // shape c
      assert( c > 0. );
   
      return pow( _u(), 1. / c );
   }
   
   double rayleigh( double a, double b ) { // Rayleigh
                                           // location a, scale b
      assert( b > 0. );
   
      return a + b * sqrt( -log( _u() ) );
   }
   
   double studentT( int df ) { // Student's T
                               // degres of freedom df
      assert( df >= 1 );
   
      return normal() / sqrt( chiSquare( df ) / df );
   }
   
   double triangular( double xMin = 0.,     // Triangular
                      double xMax = 1.,     // with default interval [0,1)
                      double c    = 0.5 ) { // and default mode 0.5

      assert( xMin < xMax && xMin <= c && c <= xMax );
      
      double p = _u(), q = 1. - p;
      
      if ( p <= ( c - xMin ) / ( xMax - xMin ) )
         return xMin + sqrt( ( xMax - xMin ) * ( c - xMin ) * p );
      else
         return xMax - sqrt( ( xMax - xMin ) * ( xMax - c ) * q );
   }
   
   double uniform( double xMin = 0., double xMax = 1. ) { // Uniform
                                                          // on [xMin,xMax)
      assert( xMin < xMax );
   
      return xMin + ( xMax - xMin ) * _u();
   }
   
   double userSpecified(                // User-Specified Distribution
        double( *usf )(                 // pointer to user-specified function
             double,                    // x
             double,                    // xMin
             double ),                  // xMax
        double xMin, double xMax,       // function domain
        double yMin, double yMax ) {    // function range

      assert( xMin < xMax && yMin < yMax );
   
      double x, y, areaMax = ( xMax - xMin ) * ( yMax - yMin ); 

      // acceptance-rejection method
   
      do {   
         x = uniform( 0.0, areaMax ) / ( yMax - yMin ) + xMin;  
         y = uniform( yMin, yMax );
   
      } while ( y > usf( x, xMin, xMax ) );
   
      return x;
   }
   
   double weibull( double a, double b, double c ) { // Weibull
                                                    // location a, scale b,
      assert( b > 0. && c > 0. );                   // shape c
   
      return a + b * pow( -log( _u() ), 1. / c );
   }
                   
// Discrete Distributions

   bool bernoulli( double p = 0.5 ) { // Bernoulli Trial
   
      assert( 0. <= p && p <= 1. );
   
      return _u() < p;
   }
   
   int binomial( int n, double p ) { // Binomial

      assert( n >= 1 && 0. <= p && p <= 1. );
   
      int sum = 0;
      for ( int i = 0; i < n; i++ ) sum += bernoulli( p );
      return sum;
   }
   
   int geometric( double p ) { // Geometric

      assert( 0. < p && p < 1. );

      return int( log( _u() ) / log( 1. - p ) );
   }
   
   int hypergeometric( int n, int N, int K ) {          // Hypergeometric
                                                        // trials n, size N,
      assert( 0 <= n && n <= N && N >= 1 && K >= 0 );   // successes K
      
      int count = 0;
      for ( int i = 0; i < n; i++, N-- ) {
   
         double p = double( K ) / double( N );
         if ( bernoulli( p ) ) { count++; K--; }
      }
      return count;
   }
   
   void multinomial( int    n,            // Multinomial
                     double p[],          // trials n, probability vector p,
                     int    count[],      // success vector count,
                     int    m ) {         // number of disjoint events m

      assert( m >= 2 );   // at least 2 events
      double sum = 0.;
      for ( int bin = 0; bin < m; bin++ ) sum += p[ bin ];    // probabilities
      assert( sum == 1. );                                    // must sum to 1
      
      for ( int bin = 0; bin < m; bin++ ) count[ bin ] = 0;   // initialize
   
      // generate n uniform variates in the interval [0,1) and bin the results
   
      for ( int i = 0; i < n; i++ ) {

         double lower = 0., upper = 0., u = _u();

         for ( int bin = 0; bin < m; bin++ ) {

         // locate subinterval, which is of length p[ bin ],
         // that contains the variate and increment the corresponding counter
      
            lower = upper;
            upper += p[ bin ];
            if ( lower <= u && u < upper ) { count[ bin ]++; break; }
         }
      }
   }
   
   int negativeBinomial( int s, double p ) { // Negative Binomial
                                             // successes s, probability p
      assert( s >= 1 && 0. < p && p < 1. );
   
      int sum = 0;
      for ( int i = 0; i < s; i++ ) sum += geometric( p );
      return sum;
   }
   
   int pascal( int s, double p ) { // Pascal
                                   // successes s, probability p
      return negativeBinomial( s, p ) + s;
   }
   
   int poisson( double mu ) { // Poisson

      assert ( mu > 0. );
   
      double a = exp( -mu );
      double b = 1.;
   
      int i;
      for ( i = 0; b >= a; i++ ) b *= _u();   
      return i - 1;
   }
   
   int uniformDiscrete( int i, int j ) { // Uniform Discrete
                                         // inclusive i to j
      assert( i < j );

      return i + int( ( j - i + 1 ) * _u() );
   }

// Empirical and Data-Driven Distributions

   double empirical( void ) { // Empirical Continuous

      static std::vector< double > x, cdf;
      static int              n;
      static bool             init = false;
   
      if ( !init ) {
         std::ifstream in( "empiricalDistribution" );
         if ( !in ) {
            std::cerr << "Cannot open \"empiricalDistribution\" input file" << std::endl;
            exit( 1 );
         }
         double value, prob;
         while ( in >> value >> prob ) {   // read in empirical distribution
            x.push_back( value );
            cdf.push_back( prob );
         }
         n = x.size();
         init = true;
      
         // check that this is indeed a cumulative distribution

         assert( 0. == cdf[ 0 ] && cdf[ n - 1 ] == 1. );
         for ( int i = 1; i < n; i++ ) assert( cdf[ i - 1 ] < cdf[ i ] );
      }

      double p = _u();
      for ( int i = 0; i < n - 1; i++ )
         if ( cdf[ i ] <= p && p < cdf[ i + 1 ] )
            return x[ i ] + ( x[ i + 1 ] - x[ i ] ) *
                            ( p - cdf[ i ] ) / ( cdf[ i + 1 ] - cdf[ i ] );
      return x[ n - 1 ];
   }
   
   int empiricalDiscrete( void ) { // Empirical Discrete

      static std::vector< int >    k;
      static std::vector< double > f[ 2 ];   // f[ 0 ] is pdf and f[ 1 ] is cdf
      static double           max;
      static int              n;
      static bool             init = false;
   
      if ( !init ) {
         std::ifstream in ( "empiricalDiscrete" );
         if ( !in ) {
            std::cerr << "Cannot open \"empiricalDiscrete\" input file" << std::endl;
            exit( 1 );
         }
         int value;
         double freq;
         while ( in >> value >> freq ) {   // read in empirical data
            k.push_back( value );
            f[ 0 ].push_back( freq );
         }
         n = k.size();
         init = true;
   
         // form the cumulative distribution

         f[ 1 ].push_back( f[ 0 ][ 0 ] );
         for ( int i = 1; i < n; i++ )
            f[ 1 ].push_back( f[ 1 ][ i - 1 ] + f[ 0 ][ i ] );

         // check that the integer points are in ascending order

         for ( int i = 1; i < n; i++ ) assert( k[ i - 1 ] < k[ i ] );
      
         max = f[ 1 ][ n - 1 ];
      }
   
      // select a uniform variate between 0 and the max value of the cdf

      double p = uniform( 0., max );

      // locate and return the corresponding index

      for ( int i = 0; i < n; i++ ) if ( p <= f[ 1 ][ i ] ) return k[ i ];
      return k[ n - 1 ];
   }
   
   double sample( bool replace = true ) { // Sample w or w/o replacement from a
                                          // distribution of 1-D data in a file
      static std::vector< double > v;     // vector for sampling with replacement
      static bool init = false;           // flag that file has been read in
      static int n;                       // number of data elements in the file
      static int index = 0;               // subscript in the sequential order
   
      if ( !init ) {
         std::ifstream in( "sampleData" );
         if ( !in ) {
            std::cerr << "Cannot open \"sampleData\" file" << std::endl;
            exit( 1 );
         }
         double d;
         while ( in >> d ) v.push_back( d );
         in.close();
         n = v.size();
         init = true;
         if ( replace == false ) {   // sample without replacement
         
            // shuffle contents of v once and for all
            // Ref: Knuth, D. E., The Art of Computer Programming, Vol. 2:
            //      Seminumerical Algorithms. London: Addison-Wesley, 1969.

            for ( int i = n - 1; i > 0; i-- ) {
               int j = int( ( i + 1 ) * _u() );
               std::swap( v[ i ], v[ j ] );
            }
         }
      }

      // return a random sample

      if ( replace )                                // sample w/ replacement
         return v[ uniformDiscrete( 0, n - 1 ) ];
      else {                                        // sample w/o replacement
         assert( index < n );                       // retrieve elements
         return v[ index++ ];                       // in sequential order
      }
   }
   
   void sample( double x[], int ndim ) { // Sample from a given distribution
                                         // of multi-dimensional data
      static const int N_DIM = 6;
      assert( ndim <= N_DIM );
   
      static std::vector< double > v[ N_DIM ];
      static bool init = false;
      static int n;
   
      if ( !init )  {
         std::ifstream in( "sampleData" );
         if ( !in ) {
            std::cerr << "Cannot open \"sampleData\" file" << std::endl;
            exit( 1 );
         }
         double d;
         while ( !in.eof() ) {
            for ( int i = 0; i < ndim; i++ ) {
               in >> d;
               v[ i ].push_back( d );
            }
         }
         in.close();
         n = v[ 0 ].size();
         init = true;
      }
      int index = uniformDiscrete( 0, n - 1 );
      for ( int i = 0; i < ndim; i++ ) x[ i ] = v[ i ][ index ];
   }

   // comparison functor for determining the neighborhood of a data point

   struct dSquared :
   public std::binary_function< point2d, point2d, bool > {
      bool operator()( point2d p, point2d q ) {
         return p.x * p.x + p.y * p.y < q.x * q.x + q.y * q.y;
      }
   };

   point2d stochasticInterpolation( void ) { // Stochastic Interpolation

   // Refs: Taylor, M. S. and J. R. Thompson, Computational Statistics & Data 
   //       Analysis, Vol. 4, pp. 93-101, 1986; Thompson, J. R., Empirical Model
   //       Building, pp. 108-114, Wiley, 1989; Bodt, B. A. and M. S. Taylor,
   //       A Data Based Random Number Generator for A Multivariate Distribution
   //       - A User's Manual, ARBRL-TR-02439, BRL, APG, MD, Nov. 1982.

      static std::vector<point2d> data;
      static point2d              min, max;
      static int                  m;
      static double               lower, upper;
      static bool                 init = false;

      if ( !init ) {
         std::ifstream in( "stochasticData" );
         if ( !in ) {
            std::cerr << "Cannot open \"stochasticData\" input file" << std::endl;
            exit( 1 );
         }
   
         // read in the data and set min and max values
   
         min.x = min.y = FLT_MAX;
         max.x = max.y = FLT_MIN;
         point2d p;
         while ( in >> p.x >> p.y ) {

            min.x = ( p.x < min.x ? p.x : min.x );
            min.y = ( p.y < min.y ? p.y : min.y );
            max.x = ( p.x > max.x ? p.x : max.x );
            max.y = ( p.y > max.y ? p.y : max.y );
      
            data.push_back( p );
         }
         in.close();
         init = true;
   
         // scale the data so that each dimension will have equal weight

         for ( unsigned int i = 0; i < data.size(); i++ ) {
       
            data[ i ].x = ( data[ i ].x - min.x ) / ( max.x - min.x );
            data[ i ].y = ( data[ i ].y - min.y ) / ( max.y - min.y );
         }
   
         // set m, the number of points in a neighborhood of a given point
   
         m = data.size() / 20;       // 5% of all the data points
         if ( m < 5  ) m = 5;        // but no less than 5
         if ( m > 20 ) m = 20;       // and no more than 20
   
         lower = ( 1. - sqrt( 3. * ( double( m ) - 1. ) ) ) / double( m );
         upper = ( 1. + sqrt( 3. * ( double( m ) - 1. ) ) ) / double( m );
      }
   
      // uniform random selection of a data point (with replacement)
      
      point2d origin = data[ uniformDiscrete( 0, data.size() - 1 ) ];

      // make this point the origin of the coordinate system

      for ( unsigned int n = 0; n < data.size(); n++ ) data[ n ] -= origin;
      
      // sort the data with respect to its distance (squared) from this origin
      
      std::sort( data.begin(), data.end(), dSquared() );
      
      // find the mean value of the data in the neighborhood about this point
      
      point2d mean;
      mean.x = mean.y = 0.;
      for ( int n = 0; n < m; n++ ) mean += data[ n ];
      mean /= double( m );
   
      // select a random linear combination of the points in this neighborhood

      point2d p;
      p.x = p.y = 0.;
      for ( int n = 0; n < m; n++ ) {
         
         double rn;
         if ( m == 1 ) rn = 1.;
         else          rn = uniform( lower, upper );
         p.x += rn * ( data[ n ].x - mean.x );
         p.y += rn * ( data[ n ].y - mean.y );
      }
      
      // restore the data to its original form

      for ( unsigned int n = 0; n < data.size(); n++ ) data[ n ] += origin;
      
      // use mean and original point to translate the randomly-chosen point

      p += mean;
      p += origin;

      // scale randomly-chosen point to the dimensions of the original data
      
      p.x = p.x * ( max.x - min.x ) + min.x;
      p.y = p.y * ( max.y - min.y ) + min.y;

      return p;
   }

// Multivariate Distributions

   point2d bivariateNormal( double muX    = 0.,   // Bivariate Gaussian
                            double sigmaX = 1.,
                            double muY    = 0., 
                            double sigmaY = 1. ) {

      assert( sigmaX > 0. && sigmaY > 0. );
   
      point2d p;
      p.x = normal( muX, sigmaX );
      p.y = normal( muY, sigmaY );
      return p;
   }
   
   point2d bivariateUniform( double xMin = -1.,    // Bivariate Uniform
                             double xMax =  1.,
                             double yMin = -1.,
                             double yMax =  1. ) {

      assert( xMin < xMax && yMin < yMax );
   
      double x0 = 0.5 * ( xMin + xMax );
      double y0 = 0.5 * ( yMin + yMax );
      double a  = 0.5 * ( xMax - xMin );
      double b  = 0.5 * ( yMax - yMin );
      double x, y;
   
      do {
         x = uniform( -1., 1. );
         y = uniform( -1., 1. );
      
      } while( x * x + y * y > 1. );
      
      point2d p;
      p.x = x0 + a * x;
      p.y = y0 + b * y;
      return p;
   }
   
   point2d corrNormal( double r,              // Correlated Normal
                       double muX    = 0.,
                       double sigmaX = 1.,
                       double muY    = 0.,
                       double sigmaY = 1. ) {

      assert( -1. <= r && r <= 1. &&          // bounds on correlation coeff
              sigmaX > 0. && sigmaY > 0. );   // positive std dev
   
      double x = normal();
      double y = normal();
   
      y = r * x + sqrt( 1. - r * r ) * y;     // correlate the variables
   
      point2d p;
      p.x = muX + sigmaX * x;                 // translate and scale
      p.y = muY + sigmaY * y;
      return p;
   }
   
   point2d corrUniform( double r,        // Correlated Uniform
                        double xMin = 0.,
                        double xMax = 1.,
                        double yMin = 0.,
                        double yMax = 1. ) {

      assert( -1. <= r && r <= 1. &&          // bounds on correlation coeff
              xMin < xMax && yMin < yMax );   // bounds on domain

      double x0 = 0.5 * ( xMin + xMax );
      double y0 = 0.5 * ( yMin + yMax );
      double a  = 0.5 * ( xMax - xMin );
      double b  = 0.5 * ( yMax - yMin );
      double x, y;
   
      do {
         x = uniform( -1., 1. );
         y = uniform( -1., 1. );
      
      } while ( x * x + y * y > 1. );
   
      y = r * x + sqrt( 1. - r * r ) * y;   // correlate the variables
   
      point2d p;
      p.x = x0 + a * x;                     // translate and scale
      p.y = y0 + b * y;
      return p;
   }

private:

   // function returns an associative array where each element of the vector is the key and the rank is the value
   inline std::map<double,int> _rank( std::vector<double> v ) {   // NB: pass a copy of the vector, not a reference

      std::sort( v.begin(), v.end() );   // sort the vector in ascending order
      std::map<double,int> r;
      for ( unsigned int i = 0; i < v.size(); i++ ) r[v[i]] = i;   // map the element to its rank
      return r;
   }

public:

   // correlate two distributions without changing the marginal distributions
   // (Thanks to Dr. Joseph Collins for describing this technique.  For an understanding of the theory,
   // and more general cases of any number of distributions, see ARL-TR-XXX and references cited therein.)
   void corrDist( std::vector<double>& dist1, std::vector<double>& dist2, double rankCorr ) {   // the two distributions are reordered to induce the correlation

      std::vector<double> t1( dist1 ), t2( dist2 );   // copy the two distributions
      std::sort( t1.begin(), t1.end() );   // sort the copies in ascending order
      std::sort( t2.begin(), t2.end() );
   
      double x, y, c = 2. * sin( rankCorr * M_PI / 6. ), s = sqrt( 1. - c * c );
      const int N = dist1.size();
      std::vector<double> u(N), v(N);
      for ( int i = 0; i < N; i++ ) {   // generate two correlated vectors with the given rank correlation

         x = normal();
         y = normal();
         y = c * x + s * y;   // perform a rotation to induce the correlation
         u[i] = x ;
         v[i] = y;
      }
      std::map<double,int> rank_u = _rank( u );   // generate maps from the values to the corresponding ranks
      std::map<double,int> rank_v = _rank( v );
   
      for ( int i = 0; i < N; i++ ) {   // apply these maps as index functions to the sorted distributions

         dist1[i] = t1[ rank_u[ u[i] ] ];
         dist2[i] = t2[ rank_v[ v[i] ] ];
      }
   }

   point3d spherical( double thMin = 0.,     // Uniform Spherical
                             double thMax = M_PI,
                             double phMin = 0.,
                             double phMax = 2. * M_PI ) {

      assert( 0. <= thMin && thMin < thMax && thMax <= M_PI &&
              0. <= phMin && phMin < phMax && phMax <= 2. * M_PI );
   
      point3d p;                       // azimuth
      p.theta = acos( uniform( cos( thMax ), cos( thMin ) ) );   // polar angle
      p.phi   = uniform( phMin, phMax );                         // azimuth
      return p;
   }
   
   void sphericalND( double x[], int n ) { // Uniform over the surface of
                                           // an n-dimensional unit sphere

      // generate a point inside the unit n-sphere by normal polar method

      double r2 = 0.;
      for ( int i = 0; i < n; i++ ) {
         x[ i ] = normal();
         r2 += x[ i ] * x[ i ];
      }
   
      // project the point onto the surface of the unit n-sphere by scaling

      const double A = 1. / sqrt( r2 );
      for ( int i = 0; i < n; i++ ) x[ i ] *= A;
   }

// Number Theoretic Distributions

   double avoidance( void ) { // Maximal Avoidance (1-D)
                              // overloaded for convenience
      double x[ 1 ];
      avoidance( x, 1 );
      return x[ 0 ];
   }
   
   void avoidance( double x[], unsigned int ndim ) { // Maximal Avoidance (N-D)

      static const unsigned int MAXBIT = 30;
      static const unsigned int MAXDIM = 6;
   
      assert( ndim <= MAXDIM );

      static unsigned long ix[ MAXDIM + 1 ] = { 0 };
      static unsigned long *u[ MAXBIT + 1 ];
      static unsigned long mdeg[ MAXDIM + 1 ] = { // degree of
         0,                                       // primitive polynomial
         1, 2, 3, 3, 4, 4
      };
      static unsigned long p[ MAXDIM + 1 ] = {   // decimal encoded
         0,                                      // interior bits
         0, 1, 1, 2, 1, 4
      };
      static unsigned long v[ MAXDIM * MAXBIT + 1 ] = {
          0,
          1,  1, 1,  1,  1,  1,
          3,  1, 3,  3,  1,  1,
          5,  7, 7,  3,  3,  5,
         15, 11, 5, 15, 13,  9
      };

      static double fac;
      static int in = -1;
      unsigned int j, k;
      unsigned long i, m, pp;
      
      if ( in == -1 ) {
         in = 0;
         fac = 1. / ( 1L << MAXBIT );
         for ( j = 1, k = 0; j <= MAXBIT; j++, k += MAXDIM ) u[ j ] = &v[ k ];
         for ( k = 1; k <= MAXDIM; k++ ) {
            for ( j = 1; j <= mdeg[ k ]; j++ ) u[ j ][ k ] <<= ( MAXBIT - j );
            for ( j = mdeg[ k ] + 1; j <= MAXBIT; j++ ) {
               pp = p[ k ];
               i = u[ j - mdeg[ k ] ][ k ];
               i ^= ( i >> mdeg[ k ] );
               for ( int n = mdeg[ k ] - 1; n >= 1; n-- ) {
                  if ( pp & 1 ) i ^= u[ j - n ][ k ];
                  pp >>= 1;
               }
               u[ j ][ k ] = i;
            }
         }
      }
      m = in++;
      for ( j = 0; j < MAXBIT; j++, m >>= 1 ) if ( !( m & 1 ) ) break;
      if ( j >= MAXBIT ) exit( 1 );
      m = j * MAXDIM;
      for ( k = 0; k < ndim; k++ ) {
         ix[ k + 1 ] ^= v[ m + k + 1 ];
         x[ k ] = ix[ k + 1 ] * fac;
      }
   }
                             
private:

   Generator<Typename> *_gen;
   static const unsigned int N_BITS = CHAR_BIT * sizeof( Typename );   // number of bits (32 or 64)
   
   long double _u( void ) {
   
      if ( N_BITS == 32 ) return _gen->rng32_01();
      else                return _gen->rng64_01();
   }
   
   static double _parabola( double x, double xMin, double xMax ) { // parabola

      if ( x < xMin || x > xMax ) return 0.0;
   
      double a    = 0.5 * ( xMin + xMax );   // location parameter
      double b    = 0.5 * ( xMax - xMin );   // scale parameter
      double yMax = 0.75 / b;
      
      return yMax * ( 1. - ( x - a ) * ( x - a ) / ( b * b ) );
   }
}; // Random class
} // rnd namespace

#endif // RANDOM_H
