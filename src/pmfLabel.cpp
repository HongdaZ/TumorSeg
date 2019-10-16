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


        return R_NilValue;
    }
} // extern "C"
