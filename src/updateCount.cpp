#include <R.h>
#include <Rinternals.h>
#include <omp.h>

extern "C" {
    SEXP updateCount( SEXP count, SEXP pred ) {
        double *ptr_count = REAL( count );
        double *ptr_pred = REAL( pred );
        int n = length( pred );
        int label = 0;
        for( int i = 0; i < n; ++ i ) {
            label = ptr_pred[ i ];
            ptr_count[ 4 * i + label - 1 ] += 1;
        }
        return R_NilValue;
    }

} // extern "C"
