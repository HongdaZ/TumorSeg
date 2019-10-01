#include <R.h>
#include <Rinternals.h>
#include <omp.h>

extern "C" {
    SEXP logZ( SEXP patient, SEXP alpha,  SEXP beta, SEXP gamma, SEXP s,
               SEXP p ) {

        return R_NilValue;
    }
} // extern "C"
