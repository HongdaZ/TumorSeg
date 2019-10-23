#include <R.h>
#include <Rinternals.h>
#include <omp.h>
#include "helper.h"
extern "C" {

    SEXP logZS( SEXP new_patient, SEXP pw, SEXP hoi, SEXP alpha_star,
                SEXP beta_star, SEXP gamma_star,
                SEXP beta_sum, SEXP gamma_sum,
                SEXP p  ) {
        double *alpha = REAL( alpha_star );
        double *beta = REAL( beta_star );
        double gamma = REAL( gamma_star )[ 0 ];
        double *b_sum = REAL( beta_sum );
        double *g_sum = REAL( gamma_sum );
        double res_b[ 3 ] = { 0 };
        double res_g = 0;

        int n_thread = REAL( p )[ 0 ];
        double res = 0;
        double *local_pw = REAL( pw );
        double *local_hoi = REAL( hoi );

        for( int i = 0; i < 8; ++ i ) {
            SEXP seq = VECTOR_ELT( new_patient, i );
            SEXP elmt = getListElement( seq, "neighbor_label" );
            SEXP pred_seg = getListElement( seq, "pred_seg" );
            double *ptr_nb_label = REAL( elmt );
            double *ptr_prd = REAL( pred_seg );
            SEXP dim = getAttrib( elmt, R_DimSymbol );
            int nc = INTEGER( dim )[ 1 ];
            double sum = 0;
            double sum_b[ 3 ] = { 0 };
            double sum_g = 0;

#pragma omp parallel for num_threads( n_thread ) default( none )          \
            firstprivate( alpha, beta, gamma, local_pw,                   \
                          local_hoi, nc, ptr_nb_label, ptr_prd )          \
                reduction( +: sum, sum_b[ : 3 ], sum_g )
            for( int j = 0; j < nc; ++ j ) {
                double sum_exp = 0;
                for( int k = 0; k < 4; ++ k ) {
                    int count1[ 3 ] = { 0 }; // count for beta's
                    int count2[ 8 ][ 5 ] = { 0 }; // count for labels (0 ~ 4)
                    int count3[ 8 ] = { 0 }; // count for gamma
                    for( int l = 0; l < 26; ++ l ) {
                        if( ptr_nb_label[ 26 * j + l ] != ( k + 1 ) &
                            ptr_nb_label[ 26 * j + l ] != 0 ) {
                            int dist = local_pw[ l ];
                            ++ count1[ dist - 1 ];
                        }
                    }
                    for( int m = 0; m < 8; ++ m ) {
                        ++ count2[ m ][ k + 1 ];
                        for( int o = 0; o < 7; ++ o ) {
                            int index = local_hoi[ 7 * m + o ];
                            int label = ptr_nb_label[ 26 * j + index - 1 ];
                            // label = 0, 1, 2, 3, 4
                            ++ count2[ m ][ label ];
                        }
                        for( int o = 0; o < 4; ++ o ) {
                            if( count2[ m ][ o + 1 ] > 0 ) {
                                ++ count3[ m ];
                            }
                        }
                    }
                    double log_p = 0;
                    log_p -= alpha[ k ];
                    for( int l = 0; l < 3; ++ l ) {
                        log_p -= beta[ l ] * count1[ l ];
                        if( ptr_prd[ j ] == ( k + 1 ) ) {
                            sum_b[ l ] += count1[ l ];
                        }

                    }
                    for( int l = 0; l < 8; ++ l ) {
                        log_p -= gamma * log( count3[ l ] );
                        if( ptr_prd[ j ] == ( k + 1 ) ) {
                            sum_g += log( count3[ l ] );
                        }
                    }
                    sum_exp += exp( log_p );
                }
                sum += - log( sum_exp );
            }
            for( int j = 0; j < 3; ++ j ) {
                b_sum[ j ] += sum_b[ j ];
            }
            g_sum[ 0 ] += sum_g;
            res += sum;
        }
        SEXP result = PROTECT( allocVector( REALSXP, 1 ) );
        double *res_ptr = REAL( result );
        res_ptr[ 0 ] = res;
        UNPROTECT( 1 );
        return result;
    }
} // extern "C"
