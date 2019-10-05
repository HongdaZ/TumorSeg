#include <R.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <omp.h>

extern "C" {
    SEXP recoverInfo( SEXP A, SEXP P ) {
        int nra, nca;
        SEXP Adim, B;

        Adim = getAttrib( A, R_DimSymbol );
        nra = INTEGER( Adim )[ 0 ];
        nca = INTEGER( Adim )[ 1 ];
        B = PROTECT( allocMatrix( INTSXP, nra * 4, nca ) );

        int *mat_A = INTEGER( A );
        int *mat_B = INTEGER( B );
        double *p = REAL( P );
        int n = p[ 0 ];
        div_t divresult;
        int rem;
#pragma omp parallel for num_threads( n ) default( none ) \
        private( divresult, rem )                   \
        firstprivate( nra, nca, mat_A, mat_B )
            for( int i = 0; i < nca; ++ i ) {
                for( int j = 0; j < nra; ++ j ) {
                    int index = nra * i + j;
                    divresult = div( mat_A[ index ], pow( 10, 8 ) );
                    // 4 for betas, and gamma
                    mat_B[ 4 * index + 1 ] = divresult.quot;
                    rem = divresult.rem;
                    divresult = div( rem, pow( 10, 7 ) );
                    mat_B[ 4 * index ] = divresult.quot;
                    rem = divresult.rem;
                    divresult = div( rem, pow( 10, 6 ) );
                    mat_B[ 4 * index + 2 ] = divresult.quot;
                    mat_B[ 4 * index + 3 ] = 1;

                    for( int k = 0; k < 8; ++ k ) {
                        rem = divresult.rem;
                        divresult = div( rem, pow( 5 , 7 - k ) );
                        mat_B[ 4 * index + 3 ] *= divresult.quot;
                    }
                }
            }
        UNPROTECT( 1 );
        return B;
    }
} // extern "C"
