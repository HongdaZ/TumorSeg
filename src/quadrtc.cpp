#include <R.h>
#include <Rinternals.h>
extern "C" {
    // inner product of two vectors
    double inner( double *a, double *b ) {
        double res = 0;
        for( int i = 0; i != 4; ++i ) {
            res += a[ i ] * b[ i ];
        }
        return res;
    }
    // Quadratic form of a 4 x 4 matrix
    double quadrtc( double *mat, double *vec ) {
        double res;
        double tmp[ 4 ] = { 0, 0, 0, 0 };
        for( int j = 0; j != 4; ++ j ) {
            for( int i = 0; i != 4; ++i ) {
                tmp[ j ] += vec[ i ] * mat[ 4 * j + i ];
            }
        }
        res = inner( tmp, vec );
        return res;
    }
} // extern "C"
