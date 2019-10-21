#include <R.h>
#include <Rinternals.h>
#include <omp.h>

extern "C" {
    SEXP updateCube( SEXP cube, SEXP index, SEXP pred ) {
        double *ptr_cube = REAL( cube );
        double *ptr_index = REAL( index );
        double *ptr_pred = REAL( pred );
        int n = length( pred );
        int vec_index = 0;
        for( int i = 0; i < n; ++ i ) {
            vec_index = ptr_index[ i ];
            ptr_cube[ vec_index - 1 ] = ptr_pred[ i ];
        }
        return R_NilValue;
    }

} // extern "C"
