#include <R.h>
#include <Rinternals.h>
#include <omp.h>

extern "C" {
    SEXP logZ( SEXP patient, SEXP alpha,  SEXP beta, SEXP gamma, SEXP p ) {
        double *vec_alpha = REAL( alpha );
        double *vec_beta = REAL( beta );
        double *ptr_gamma = REAL( gamma );

        SEXP elmt = VECTOR_ELT( patient, 3 );
        int *nb_info = INTEGER( elmt );
        SEXP dim = getAttrib( elmt, R_DimSymbol );
        int nc = INTEGER( dim )[ 1 ];
        int n_thread = REAL( p )[ 0 ];
        double sum = 0;
#pragma omp parallel for num_threads( n_thread ) default( none )    \
        firstprivate( nc, nb_info, vec_alpha, vec_beta, ptr_gamma ) \
        reduction( +: sum )
        for( int i = 0; i < nc; ++ i ) {
            double exp_sum = 0;
            for( int j = 0; j < 4; ++ j ) { // 4 labels
                double inner = 0;
                inner = -vec_alpha[ j ];
                for( int k = 0; k < 3; ++ k ) {
                    inner -= vec_beta[ k ] *
                        nb_info[ 16 * i + 4 * j + k ];
                }
                inner -= ptr_gamma[ 0 ] *
                    log( nb_info[ 16 * i + 4 * j + 3 ] );
                exp_sum += exp( inner );
            }
            sum += log( exp_sum );
        }
        SEXP res = PROTECT( allocVector( REALSXP, 1 ) );
        double *res_ptr = REAL( res );
        res_ptr[ 0 ] = -sum;
        UNPROTECT( 1 );
        return res;
    }
} // extern "C"
