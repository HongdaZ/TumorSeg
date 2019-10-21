#include <R.h>
#include <Rinternals.h>
#include <omp.h>

extern "C" {
    SEXP updatePred( SEXP seg, SEXP pred ) {
        double *ptr_seg = REAL( seg );
        double *ptr_pred = REAL( pred );
        int n = length( pred );
        for( int i = 0; i < n; ++ i ) {
            ptr_seg[ i ] = ptr_pred[ i ];
        }
        return R_NilValue;
    }

} // extern "C"
