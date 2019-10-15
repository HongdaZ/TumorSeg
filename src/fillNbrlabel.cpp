#include <R.h>
#include <Rinternals.h>
#include <omp.h>

extern "C" {
    SEXP fillNbrlabel( SEXP cube, SEXP neighbor_index,  SEXP neighbor_label,
                       SEXP p ) {
        SEXP dim = getAttrib( neighbor_index, R_DimSymbol );
        double *ptr_cube = REAL( cube );
        double *ptr_index = REAL( neighbor_index );
        double *ptr_label = REAL( neighbor_label );
        int nr = INTEGER( dim )[ 0 ];
        int nc = INTEGER( dim )[ 1 ];
        int n_thread = REAL( p )[ 0 ];
#pragma omp parallel for num_threads( n_thread ) default( none ) \
        firstprivate( nr, nc, ptr_cube, ptr_index, ptr_label )
        for( int j = 0; j < nc; ++j ) {
            for( int i = 0; i < nr; ++i ) {
                int index = ptr_index[ i + nr * j ];
                if( index == 0 ) {
                    ptr_label[ i + nr * j ] = 0;
                } else {
                    ptr_label[ i + nr * j ] =
                        ptr_cube[ index - 1 ];
                }
            }
        }
        return R_NilValue;
    }
} // extern "C"
