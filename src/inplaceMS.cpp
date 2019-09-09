#include <R.h>
#include <Rinternals.h>
extern "C" {
    void inplaceMS( SEXP A, SEXP B ) {
        int nra, nca;
        SEXP Adim;
        Adim = getAttrib( A, R_DimSymbol );
        nra = INTEGER( Adim )[ 0 ];
        nca = INTEGER( Adim )[ 1 ];
        int i = 0, j = 0;
        double *L = REAL( A ), *R = REAL( B );
        for( i = 0; i != nra; ++ i ) {
            for( j = 0; j != nca; ++ j  ) {
                L[ j * nra + i ] += R[ j * nra + i ];
            }
        }

    }
} // extern "C"
