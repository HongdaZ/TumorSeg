#include <R.h>
#include <Rinternals.h>
#include <omp.h>

extern "C" {
    SEXP vecIndex( SEXP mat, SEXP dim, SEXP p ) {
        SEXP mat_dim = getAttrib( mat, R_DimSymbol );
        SEXP vec_index;
        int nc = INTEGER( mat_dim )[ 1 ];
        vec_index = PROTECT( allocVector( REALSXP, nc ) );
        double *ptr_index = REAL( vec_index );
        double *ptr_mat = REAL( mat );
        int *ptr_dim = INTEGER( dim );
        int n_thread = REAL( p )[ 0 ];
        int local_dim[ 3 ];
        for( int i = 0; i < 3; ++ i ) {
            local_dim[ i ] = ptr_dim[ i ];
        }

#pragma omp parallel for num_threads( n_thread ) default( none )          \
        firstprivate( ptr_mat, ptr_index, local_dim, nc )
        for( int i = 0; i < nc; ++ i ) {
            ptr_index[ i ] = local_dim[ 0 ] * ( ptr_mat[ i * 3 + 1 ] - 1 ) +
                ptr_mat[ i * 3 ] + ( local_dim[ 0 ] * local_dim[ 1 ] ) *
                ( ptr_mat[ i * 3 + 2 ] - 1 );
        }

        UNPROTECT( 1 );
        return vec_index;
    }
} // extern "C"
