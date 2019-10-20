#include <R.h>
#include <Rinternals.h>
#include <omp.h>

extern "C" {


    /* get the list element named str, or return NULL */

    SEXP getListElement( SEXP list, const char *str )
    {
        SEXP elmt = R_NilValue, names = getAttrib( list, R_NamesSymbol );


        for ( int i = 0; i < length(list ); i++ )
            if( strcmp( CHAR( STRING_ELT( names, i ) ), str ) == 0 ) {
                elmt = VECTOR_ELT( list, i );
                break;
            }
        return elmt;
    }
    SEXP logZS( SEXP new_patient, SEXP pw, SEXP hoi, SEXP alpha_star,
                SEXP beta_star, SEXP gamma_star, SEXP p  ) {
        double *alpha = REAL( alpha_star );
        double *beta = REAL( beta_star );
        double gamma = REAL( gamma_star )[ 0 ];
        int n_thread = REAL( p )[ 0 ];
        double res = 0;
        double *local_pw = REAL( pw );
        double *local_hoi = REAL( hoi );

        for( int i = 0; i < 8; ++ i ) {
            SEXP elmt = VECTOR_ELT( new_patient, i );
            elmt = getListElement( elmt, "neighbor_label" );
            double *ptr_nb_label = REAL( elmt );
            SEXP dim = getAttrib( elmt, R_DimSymbol );
            int nc = INTEGER( dim )[ 1 ];
            double sum = 0;

#pragma omp parallel for num_threads( n_thread ) default( none )          \
            firstprivate( alpha, beta, gamma, local_pw,                   \
                          local_hoi, nc, ptr_nb_label )                   \
                reduction( +: sum )
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
                    }
                    for( int l = 0; l < 8; ++ l ) {
                        log_p -= gamma * log( count3[ l ] );
                    }
                    sum_exp += exp( log_p );
                }
                sum += - log( sum_exp );
            }
            res += sum;
        }


        SEXP result = PROTECT( allocVector( REALSXP, 1 ) );
        double *res_ptr = REAL( result );
        res_ptr[ 0 ] = res;
        UNPROTECT( 1 );
        return result;
    }
} // extern "C"
