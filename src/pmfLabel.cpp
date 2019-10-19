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
    SEXP pmfLabel( SEXP alpha_star, SEXP beta_star, SEXP gamma_star,
                   SEXP det_sigma, SEXP inv_sigma, SEXP sigma_mu,
                   SEXP q_sigma, SEXP pw, SEXP hoi, SEXP neighbor_label,
                   SEXP sub_modality_mat,
                   SEXP prob, SEXP p ) {
        double *ptr_alpha = REAL( alpha_star );
        double *ptr_beta = REAL( beta_star );
        double *ptr_gamma = REAL( gamma_star );
        double *ptr_det = REAL( det_sigma );
        double *ptr_inv = REAL( inv_sigma );
        double *ptr_sm = REAL( sigma_mu );
        double *ptr_q = REAL( q_sigma );
        double *ptr_pw = REAL( pw );
        double *ptr_hoi = REAL( hoi );
        double *ptr_nb_label = REAL( neighbor_label );
        double *ptr_mod = REAL( sub_modality_mat );
        double *ptr_prob = REAL( prob );
        int n_thread = REAL( p )[ 0 ];

        double alpha[ 4 ];
        double beta[ 3 ];
        double gamma = ptr_gamma[ 0 ];
        double det[ 4 ];
        double inv[ 4 * 4 * 4 ];
        double sm[ 4 * 4 ];
        double q[ 4 ];
        double local_pw[ 26 ];
        double local_hoi[ 7 * 8 ];
        for( int i = 0; i < 4; ++ i ) {
            alpha[ i ] = ptr_alpha[ i ];
        }
        for( int i = 0; i < 3; ++ i ) {
            beta[ i ] = ptr_beta[ i ];
        }
        for( int i = 0; i < 4; ++ i ) {
            det[ i ] = ptr_det[ i ];
        }
        for( int i = 0; i < 64; ++ i ) {
            inv[ i ] = ptr_inv[ i ];
        }
        for( int i = 0; i < 16; ++ i ) {
            sm[ i ] = ptr_sm[ i ];
        }
        for( int i = 0; i < 4; ++ i ) {
            q[ i ] = ptr_q[ i ];
        }
        for( int i = 0; i < 26; ++ i  ) {
            local_pw[ i ] = ptr_pw[ i ];
        }
        for( int i = 0; i < 56; ++ i ) {
            local_hoi[ i ] = ptr_hoi[ i ];
        }
        SEXP dim = getAttrib( neighbor_label, R_DimSymbol );
        int nc = INTEGER( dim )[ 1 ];

#pragma omp parallel for num_threads( n_thread ) default( none )      \
        firstprivate( alpha, beta, gamma, det, inv, sm, q,            \
                      local_pw, local_hoi, ptr_nb_label, ptr_mod,     \
                      ptr_prob, nc )
        for( int i = 0; i < nc; ++ i ) {
            double sum_prob = 0;
            for( int j = 0; j < 4; ++ j ) {
                int count1[ 3 ] = { 0 }; // count for beta's
                int count2[ 8 ][ 5 ] = { 0 }; // count for labels (0 ~ 4)
                int count3[ 8 ] = { 0 }; // count for gamma
                for( int k = 0; k < 26; ++ k ) {
                    if( ptr_nb_label[ 26 * i + k ] != ( j + 1 ) &
                        ptr_nb_label[ 26 * i + k ] != 0 ) {
                        int dist = local_pw[ k ];
                        ++ count1[ dist - 1 ];
                    }
                }
                for( int m = 0; m < 8; ++ m ) {
                    ++ count2[ m ][ j + 1 ];
                    for( int o = 0; o < 7; ++ o ) {
                        int index = local_hoi[ 7 * m + o ];
                        int label = ptr_nb_label[ 26 * i + index - 1 ];
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
                log_p -= alpha[ j ];
                for( int k = 0; k < 3; ++ k ) {
                    log_p -= beta[ k ] * count1[ k ];
                }
                for( int k = 0; k < 8; ++ k ) {
                    log_p -= gamma * log( count3[ k ] );
                }
                log_p -= 1 / 2 * log( det[ j ] );
                log_p -= 1 / 2 * ( quadrtc( inv + 16 * j, ptr_mod + 4 * i  ) -
                    2 * inner( ptr_mod + 4 * i, sm + 4 * j ) +  q[ j ] );
                double p = exp( log_p );
                ptr_prob[ 4 * i + j ] = p;
                sum_prob += p;
            }
            for( int j = 0; j < 4; ++ j ) {
                ptr_prob[ 4 * i + j ] /= sum_prob;
            }
        }
        return R_NilValue;
    }
} // extern "C"
