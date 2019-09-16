#include <R.h>
#include <Rinternals.h>
extern "C" {
    SEXP changeLabel( SEXP A, SEXP B ) {
        int na = length( A ), nb = length( B );
        int *array = INTEGER( A );
        double *vec = REAL( B );
        for( int i = 0; i != na; ++ i ) {
            for( int j = 0; j != nb; ++ j ) {
                if( array[ i ] == vec[ j ] ) {
                    array[ i ] = j + 1;
                    break;
                }
            }
        }
        return R_NilValue;
    }
} // extern "C"
