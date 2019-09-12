#include <R.h>
#include <Rinternals.h>
extern "C" {
    SEXP inplaceMS( SEXP A, SEXP B ) {
        int nra, nca;
        SEXP Adim;
        Adim = getAttrib( A, R_DimSymbol );
        nra = INTEGER( Adim )[ 0 ];
        nca = INTEGER( Adim )[ 1 ];
        int i = 0, last = nra * nca;
        double *L = REAL( A ), *R = REAL( B );
        for( i = 0; i != last ; ++ i ) {

            L[ i ] += R[ i ];
        }
        return A;
    }
} // extern "C"
