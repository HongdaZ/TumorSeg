#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <omp.h>
extern "C" {
    SEXP recoverInfo( SEXP A, SEXP B, SEXP C,SEXP P ) {

        int n_a = length( A );
        int *vec_A = INTEGER( A );
        int *vec_B = INTEGER( B );
        int *vec_C = INTEGER( C );
        double *p = REAL( P );
        int n = p[ 0 ];
#pragma omp parallel for num_threads( n ) default( none ) \
        firstprivate(n_a, vec_A, vec_B, vec_C)
        for( int i = 0; i < n_a; ++ i ) {
            div_t divresult;
            for( int j = 0; j < 1000; ++ j ) {
                divresult = div( vec_A[ i ], 10 );
                vec_B[ i ] = divresult.quot;
                vec_C[ i ] = divresult.rem;
            }

        }
        return R_NilValue;
    }
} // extern "C"
